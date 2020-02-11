import base64
import os
import pickle
import random
import time
import zlib
from io import BytesIO
from math import ceil

import numpy as np
import requests
from PIL import Image

DEBUG = False


def trim_bytes(numpy_array, num_bytes):
    """Quick way of removing the smallest n bits from an array of bytes."""
    return numpy_array >> num_bytes


def split_input(value, n):
    return [value[i : i + n] for i in range(0, len(value), n)]


def encode_input(value):
    return zlib.compress(pickle.dumps(value, protocol=0))


def decode_input(value):
    val_bytes = value.encode("latin-1")
    return pickle.loads(zlib.decompress(val_bytes))


def _print(message, value=None, indent=0):
    if DEBUG:
        if value is None:
            print(f"{' ' * indent}{message}")
        else:
            print(f"{' ' * indent}{message}: {value}")


class Steganography(object):

    MAX_LEN = 256 ** 256 - 1
    _MARKER = "00101110001"  # Must be 11 bytes

    def __init__(self, input_data, original_data=None):

        self.i = encode_input(input_data)
        self.o_ignore = original_data is None
        self.o = np.asarray(original_data)

        # Detect whether the input should be compressed
        self.enc = True
        if isinstance(input_data, str):
            _print(f"Encoded length: {len(self.i)}")
            _print(f"Normal length: {len(input_data)}")
            if len(self.i) > len(input_data):
                self.i = input_data
                self.enc = False
        if isinstance(self.i[0], int):
            self.i_ord = [i for i in self.i]
        else:
            self.i_ord = [ord(i) for i in self.i]
        self.i_len = len(self.i_ord)
        if not self.o_ignore:
            self.o_len = len(self.o)

        # Basic error checking
        if self.i_len > self.MAX_LEN:
            raise OverflowError("input data is too large")
        if not self.i_len:
            raise ValueError("input data is empty")

    def generate_header(self):

        # Find how much space it is to store the length of data
        input_len_binary = "{0:b}".format(self.i_len)
        input_len_padding = len(input_len_binary) // 8
        input_len_binary = input_len_binary.zfill(8 * (input_len_padding + 1))

        marker = split_input(self._MARKER, 8)
        if len(marker[0]) != 8 or len(marker[1]) != 3:
            raise ValueError("incorrect marker size")
        header = [int(marker[0], 2), None, input_len_padding]
        header += [int(i, 2) for i in split_input(input_len_binary, 8)]

        # Find how many bytes per colour are needed
        if self.o_ignore:
            i = 7

        else:
            bytes_available = self.o_len - len(header)

            if bytes_available < self.i_len:
                i = 7
                print("Not enough space to encode data, reverting to default encoding")
                # raise ValueError('not enough space to encode data')

            else:
                for i in range(8):
                    if self.i_len * (8 - i) < bytes_available:
                        break

        bytes_per_colour = i + 1

        header[1] = int(marker[1] + str(int(self.enc)) + "{0:04b}".format(i + 1), 2)

        _print("Bytes per colour", bytes_per_colour, indent=1)
        _print("Input length", self.i_len, indent=1)
        _print("Compressed", bool(self.enc), indent=1)
        _print("Header bytes", len(header), indent=1)
        if not self.o_ignore:
            _print("Total bytes available", self.o_len, indent=1)
            _print("Remaining bytes", bytes_available, indent=1)

        return header, bytes_per_colour

    def encode(self):

        # Generate header
        _print("Generating header...")
        data, bytes_per_colour = self.generate_header()

        # Convert data
        if not self.o_ignore:
            _print("Converting input data to binary...")
            trimmed_original = trim_bytes(self.o, bytes_per_colour)
            format_binary = "{0:0" + "{}".format(8 - bytes_per_colour) + "b}"
            original_bytes = [format_binary.format(i) for i in trimmed_original]

        new_bytes = split_input(
            "".join(f"{i:08b}" for i in self.i_ord), bytes_per_colour
        )
        new_bytes[-1] += "0" * (bytes_per_colour - len(new_bytes[-1]))

        # Increase length until both inputs match
        _print("Matching length...")
        if bytes_per_colour == 8:
            required_length = len(data)
        else:
            required_length = len(original_bytes) - len(data)

        while True:
            new_byte_len = len(new_bytes)
            if new_byte_len >= required_length:
                break
            extra_length = min(new_byte_len, required_length - new_byte_len)
            new_bytes += new_bytes[:extra_length]
            _print(
                f"Increased length from {new_byte_len} to {new_byte_len + extra_length}",
                indent=1,
            )

        _print("Completed encoding")
        if bytes_per_colour == 8:
            return data + [int(i, 2) for i in new_bytes]
        else:
            return data + [
                int(a + b, 2) for a, b in zip(original_bytes[len(data) :], new_bytes)
            ]

    @classmethod
    def read_header(self, data):

        initial_data = f"{data[1]:08b}"
        marker = "{0:08b}".format(data[0]) + initial_data[0:3]
        if marker != self._MARKER:
            raise ValueError("not a valid file")
        bytes_per_colour = int(initial_data[4:8], 2)
        encoded = int(initial_data[3], 2)

        num_bytes = int("".join("{0:08b}".format(i) for i in data[3 : data[2] + 4]), 2)

        _print("Bytes per colour", bytes_per_colour, indent=1)
        _print("Input length", num_bytes, indent=1)
        _print("Compressed", bool(encoded), indent=1)

        return data[data[2] + 4 :], bytes_per_colour, data[2], num_bytes, encoded

    @classmethod
    def decode(self, data):

        _print("Reading header...")
        data, bytes_per_colour, data_len, num_bytes, encoded = self.read_header(data)

        _print("Decoding data...")
        encoded_data = split_input(
            "".join(f"{i:08b}"[8 - bytes_per_colour :] for i in data), 8
        )[:num_bytes]
        decoded_data = "".join(chr(int(i, 2)) for i in encoded_data)

        _print("Completed decoding")
        if encoded:
            return decode_input(decoded_data)
        else:
            return decoded_data


class ImageHelper(object):
    def __init__(self, path=None):
        self.path = path
        self.image = None
        self.size = None
        self.fsize = None

    def read_file(self, optional_path=None):
        _print("Reading file...")
        path = optional_path or self.path
        im = Image.open(path).convert("RGB")
        self.image = np.asarray(im).ravel()
        self.size = im.size
        _print("Image dimensions", self.size, indent=1)
        return self.image, self.size

    def read_url(self, optional_url=None):
        url = optional_url or self.path
        _print("Downloading URL contents...")
        r = requests.get(url)
        if r.status_code == 200:
            image_bytes = BytesIO(r.content)
        return self.read_file(image_bytes)

    @classmethod
    def generate_save_info(self, data_len, width, height, ratio):
        """Calculate final width and height based on current values."""
        if ratio is None:
            if width is not None and height is not None:
                ratio = [width, height]
                width = height = None
            else:
                ratio = [16, 9]
        elif isinstance(ratio, str) and ":" in ratio:
            ratio = list(map(int, ratio.split(":")))
        ratio_exp = pow(data_len * ratio[0] * ratio[1], 0.5)

        if width is None:
            if height is None:
                width = int(round(ratio_exp / ratio[1]))
                height = int(round(ratio_exp / ratio[0]))
            else:
                width = data_len / height
        else:
            if height is None:
                height = data_len / width

        width = int(ceil(max(1, min(width, data_len))))
        height = int(max(1, min(height, data_len)))

        while width * height < data_len:
            height += 1

        return width, height

    def save(self, data, width=None, height=None, ratio=None):
        _print("Saving image...")
        data_len = len(data)
        padding_required = data_len % 3
        if not isinstance(data, list):
            data = list(data)
        data += [random.randint(0, 255) for _ in range(padding_required)]

        _print("Calculating dimensions", indent=1)
        # Get the with and height of the data
        data_len = len(data)
        info = self.generate_save_info(data_len // 3, width, height, ratio)
        width, height = info

        # Copy row above when data ends
        start_of_row = data_len - data_len % width
        previous_row = start_of_row - width
        remaining_pixels = width * height * 3 - data_len
        data += data[previous_row + width - remaining_pixels : start_of_row]

        _print("Building image", indent=1)
        # Create image object
        im = Image.new("RGB", (width, height))
        px = im.load()

        # Write to image
        height_range = range(height)
        width_range = range(width)
        for y in height_range:
            for x in width_range:
                position = 3 * (x + y * width)
                px[x, y] = tuple(data[position : position + 3])

        # _print("Calculating filesize", indent=1)
        ## Get the filesize
        # temporary_image = StringIO.StringIO()
        # im.save(temporary_image, "PNG")
        # filesize = len(temporary_image.getvalue())
        # temporary_image.close()
        # _print("Filesize", filesize, indent=2)

        im.save(self.path, "PNG")
        _print("Saved image")

        self.image = data
        self.size = (width, height)


class ImgurHelper(object):
    def __init__(self):
        pass


class EncodeData(object):
    def __init__(self, path=None):
        self.original_data = []
        self.original_size = None
        self.path = path

    def set_background(self, data=None, image=None, image_url=None):
        if data is None:
            if image is None and url is not None:
                data, size = ImageHelper(image_url).read_url()
            elif image is not None:
                data, size = ImageHelper(image).read_file()
        elif isinstance(data, str):
            if any(data.startswith(i) for i in ("http", "https", "ftp")):
                data, size = ImageHelper(image_url).read_url()
        if data is None:
            size = None
        self.original_data = data
        self.original_size = size

    def set_data(self, data, is_file=False):
        if is_file:
            ext = os.path.splitext(data)[1]
            with open(data, "rb") as f:
                self.new_data = ["__file_enc__", ext, f.read()]
        else:
            self.new_data = data

    def encode(self, ratio=None, width=None, height=None, save=True):
        if self.new_data is None:
            raise ValueError("no data to encode, use 'set_data()'")
        result = {"Data": Steganography(self.new_data, self.original_data).encode()}
        if save:
            if self.original_size is not None:
                if width is None and height is None and ratio is None:
                    ratio = self.original_size
                elif width is None and ratio is None:
                    width = self.original_size[0]
                elif height is None and ratio is None:
                    height = self.original_size[1]
            ImageHelper(self.path).save(
                result["Data"], width=width, height=height, ratio=ratio
            )
            result["Path"] = self.path
        return result


def read_file(path):
    with open(path, "rb") as f:
        return f.read()


def save_file(path, data):
    with open(path, "wb") as f:
        return f.write(data)


class EncodeImage(object):
    """Groups together functions for encoding an image.

    Usage:
        encoded = EncodeImage('folder/file.mp3', True).set_background('http://website/im.jpg')
        encoded.save('folder/encodedimage')

        #Result: folder/encodedimage.png
    """

    def __init__(self, data, is_file):

        if is_file:
            ext = os.path.splitext(data)[1]
            self.data = ["__file_enc__", ext, read_file(data)]
        else:
            self.data = data

        self.clear_background()

    def set_background(self, path, is_url=None):
        if (
            is_url
            or is_url is None
            and any(path.startswith(i) for i in ("http", "https", "ftp"))
        ):
            self.image, self.size = ImageHelper(path).read_url()
        else:
            self.image, self.size = ImageHelper(path).read_file()
        return self

    def clear_background(self):
        self.image = self.size = None
        return self

    def save(self, path, ratio=None, width=None, height=None):
        if self.size is not None:
            if width is None and height is None and ratio is None:
                width, height = self.size
        if self.image is not None:
            result = {"Data": Steganography(self.data, self.image).encode()}
        else:
            result = {"Data": Steganography(self.data).encode()}

        filesize = ImageHelper(path).save(
            result["Data"], width=width, height=height, ratio=ratio
        )
        result["Path"] = path
        result["Size"] = filesize
        return result


class DecodeImage(object):
    """Groups together functions for decoding an image.
    Using a marker set in EncodeImage, the filetype is preserved.

    Usage:
        decoded = DecodeImage('folder/encodedimage.png')
        if decoded.is_file:
            decoded.save('folder/decodedfile')
        else:
            print decoded
    """

    def __init__(self, path, is_url=None):

        if (
            is_url
            or is_url is None
            and any(path.startswith(i) for i in ("http", "https", "ftp"))
        ):
            image, size = ImageHelper(path).read_url()
        else:
            image, size = ImageHelper(path).read_file()
        self.decode = Steganography.decode(image)
        self.is_file = False

        if isinstance(self.decode, list) and self.decode[0] == "__file_enc__":
            self.is_file = True
            self.ext = self.decode[1]
            self.decode = self.decode[2]

    def __str__(self):
        return str(self.decode)

    def save(self, path, force_original_extension=True):
        if not self.is_file:
            raise ValueError("encoded data was not a file")
        if force_original_extension:
            # This needs improvement, but works for now
            path += ".{}".format(self.ext)
        save_file(path, self.decode)


if __name__ == "__main__":
    DEBUG = True
