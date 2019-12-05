#!/opt/local/bin/python

import argparse
import os
import sys

import mdtraj as md
#### Requires the Glotzer gsd module ####
from gsd.fl import GSDFile as GSDFile
from gsd.hoomd import HOOMDTrajectory as HOOMDTrajectory
from numpy import arange, array, loadtxt

from cme_utils.analyze import autocorr as opva
from cme_utils.analyze import diffractometer
from cme_utils.manip.builder import hoomd_xml

#########################################



def main():
    def diffract_frame(xyz, box, working_dir=False, filetype=None):

        if filetype == "dcd":
            # since mdtraj assues are units are in angstroms, we need to convert
            xyz *= 10.0
            box *= 10.0

        D = diffractometer.diffractometer(working_dir=working_dir)
        print(working_dir)
        zoom = 4
        print("Warning zoom is at {}".format(zoom))
        if working_dir == False:
            D.set(zoom=zoom, n_views=10, length_scale=3.905)
        else:
            D.set(zoom=zoom, n_views=10, length_scale=3.905)
        D.load(xyz[0], [box[0][0], box[0][1], box[0][2]])
        D.prep_matrices()
        D.average()

    def load_traj(traj_file, top_file):
        print("Loading files, this may take awhile")
        if atom_names != None:
            atom_list = make_list(atom_names)
            t = md.load(traj_file, top=top_file)
            atom_index = t.top.select(atom_list)
            t_slice = t.atom_slice(atom_index)
            return t_slice
        if atom_ids != None:
            t = md.load(traj_file, top=top_file)
            t_slice = t.atom_slice(atom_ids)
            return t_slice
        else:
            t = md.load(traj_file, top=top_file)
            return t
        print("files loaded")

    def load_gsd(traj_file):
        f = GSDFile(name=traj_file, mode="rb")
        t = HOOMDTrajectory(f)
        return t

    def load_xml(top_file):
        class xmlHolder:
            """
            Class to transform and hold the xml data in 
            the same shape as the dcd trajectory information.
            """

            xyz = array([])
            box = array([])

        def select_types_ids(types, xyz):
            if atom_names != None:
                atom_list = make_list(atom_names)
                if args.exclude_atom_name:
                    mask = [i for i, x in enumerate(types) if x not in atom_names]
                if args.include_atom_name:
                    mask = [i for i, x in enumerate(types) if x in atom_names]
                updated_xyz = [x for i, x in enumerate(xyz) if i in mask]
                return array(updated_xyz)
            if atom_ids != None:
                updated_xyz = [x for i, x in enumerate(xyz) if i in atom_ids]
                return array(updated_xyz)

        a = hoomd_xml.hoomd_xml(top_file)
        t = [[]]
        xyz, box = a.get_xyz_box()
        if atom_names != None or atom_ids != None:
            types = a.config.find("type").text.split()
            xyz = select_types_ids(types, xyz)
        t[0] = xmlHolder()
        t[0].xyz = array([xyz])
        t[0].unitcell_lengths = array([box])
        return t

    def loadData(filetype=None, top_file=None, traj_file=None):
        """
        From the ftype argument uses if statements
        to decide which way to load the data.
        Requires:
            None
        Returns:
            t - A trajectory object
        Currently supports the following file types:
            dcd
            xml
            hoomdxml
            gsd
        """
        if filetype == None:
            filetype = input(
                "No filetype specified! What is it? (dcd, xml, hoomdxml, gsd) --> "
            )
        if filetype == "dcd":
            t = load_traj(traj_file, top_file)
        if filetype == "xml" or args.filetype == "hoomdxml":
            t = load_xml(top_file)
        if filetype == "gsd":
            t = load_gsd(traj_file)
        return t

    def generate_frames_array(
        last_frame=False,
        autocorr=False,
        allf=False,
        start_t=0,
        log_write_time=1e5,
        dcd_write_time=1e6,
    ):
        """
        Generates an array with the desired 
        frames to diffract.
            Requires: None
            Returns: Array
        Options are: last frame (last_frame; including xml only),
            autocerrelated frames (autocorr),
            and all frames (allf)
            specified with booleans.
        The values to calculate the autocorrelation of potential energy
        can also be passed. They are defaulted to:
            start_t = 0,
            log_write_time = 1e5,
            dcd_write_time = 1e6.
        Will default to final frame if no frame selection
        information is given.
        """
        if last_frame is False and autocorr is False and allf is False:
            print("No frame choice was given. Defaulting to final frame.")
            frames = array([-1])
            return frames
        if last_frame:
            frames = array([-1])
            return frames
        if autocorr is True:
            print(
                "You have specified autocorrelation. Calling autocorr with the values:"
            )
            print(
                "start_t = {}, log_write_time = {}, dcd_write_time = {}.".format(
                    start_t, log_write_time, dcd_write_time
                )
            )
            framedata = opva.autocorrelation(
                start_t=start_t,
                log_write_time=log_write_time,
                dcd_write_time=dcd_write_time,
            )
            frames = arange(
                int(framedata[0]),
                int(framedata[0] + framedata[1] * framedata[2] + 1),
                int(framedata[1]),
            )
            return frames
        if allf is True:
            frames = arange(0, len(t), 1)
            return frames

    def construct_dir_name(frame, dirname):
        """
        Generates a name for the output directory.
            Requires: int -- frame number
                dirname - str or None for directory name,
                    will default to 'diffract'.
            Returns: string -- directory name
        If given -1 and no custom name, 
        it uses ``diffract'' as the default name.
        Otherwise it will use diffract-{frame} or
        {customname}-{frame} as the name.
        """
        if frame == -1:
            if dirname != None:
                return dirname
            else:
                dirname = "./diffract".format(frame)
                return dirname
        else:
            if dirname != None:
                dirname = "./{}-frame-{}".format(dirname, frame)
                return dirname
            else:
                dirname = "./diffract-frame-{}".format(frame)
                return dirname

    def relate_typeid_to_type(typeid, types):
        """
        Uses the typeids, saved by hoomd, to relate
        the unique types to the typeids.
        Required:
            typeid: list of integers describing the type number
            types: list of unique types in the trajectory.
        Returns:
            full_list: List where typeids are converted to their
                types.
        """
        full_list = [types[i] for i in typeid]
        return full_list

    def filter_gsd(snapshot, exclude_atom_name=False, include_atom_name=False):
        """
        Filters the gsd snapshot to obtain only the atom types or 
        indices desired.
        Requires:
            snapshot - GSD snapshot
        Returns:
            xyz - array of xyz coordinates only for the desired 
                particles.
        NOTE: Supports either: 
            atom type exclusions or inclusions, 
            atom_typing or inclusion of atom_ids (not exclusion of id)
        """
        xyz = snapshot.particles.position
        if atom_names != None:
            types = relate_typeid_to_type(
                snapshot.particles.typeid, snapshot.particles.types
            )
            if args.exclude_atom_name:
                mask = array([i for i, x in enumerate(types) if x not in atom_names])
            if args.include_atom_name:
                mask = array([i for i, x in enumerate(types) if x in atom_names])
                print(mask)
            updated_xyz = xyz[mask]
            return updated_xyz
        if atom_ids != None:
            updated_xyz = [x for i, x in enumerate(xyz) if i in atom_ids]
            return array(updated_xyz)

    def iterate_through_desired_frames(t, frames, dirname, filetype):
        """
        Function calls the diffract for each frame number 
        given in the array of frames.
        Requires:
            t - trajectory object
            frames - array with the frame numbers to be diffracted.
        Returns:
            None
        """
        for frame in frames:
            dirname = construct_dir_name(frame, dirname)
            dirname = make_dir(dirname)
            os.chdir(dirname)
            if filetype == "gsd":
                snapshot = t[int(frame)]
                box = snapshot.configuration.box
                if atom_names != None or atom_ids != None:
                    xyz = filter_gsd(snapshot)
                else:
                    xyz = snapshot.particles.position
                diffract_frame(array([xyz]), array([box]), working_dir=True)
            else:
                diffract_frame(
                    t[frame].xyz,
                    t[frame].unitcell_lengths,
                    working_dir=True,
                    filetype=filetype,
                )
            os.chdir("../")

    def make_list(atom_names):
        list_of_atoms = ""
        if args.include_atom_name is True:
            for atom_name in atom_names:
                list_of_atoms += "name == '" + str(atom_name) + "' or "
            list_of_atoms = list_of_atoms[:-4]  # removes trailing 'or' and white space
            return list_of_atoms
        if args.exclude_atom_name is True:
            for atom_name in atom_names:
                list_of_atoms += "name != '" + str(atom_name) + "' and "
            list_of_atoms = list_of_atoms[:-5]  # removes trailing 'and' and white space
            return list_of_atoms

    def make_dir(directory):
        # checks to see if the directory exists first
        if os.path.exists(directory):
            directory = input("Directory name already in use. Choose another --> ")
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory

    parser = argparse.ArgumentParser(
        description="Takes in traj files to generate diffraction data in /diffract by default. The order of the flags is not important. If -all, -a or -l is not specified, it will default to -l."
    )

    ##### Arguements for which file types needed. ######
    parser.add_argument(
        "-ftype",
        "--filetype",
        help="Tells which file type to use to extract position/box information. Currently supported are:\nxml and hoomdxml --- must also use '-tp' flag to specify xml file e.g. '-ftype xml -tp restart.xml'\ndcd --- must also use '-tj' and '-tp'.\ngsd --- must also use '-tj'",
    )
    parser.add_argument(
        "-tj",
        "--traj_file",
        help="this file provides trajectory information, ie a dcd or gsd file",
    )
    parser.add_argument(
        "-tp",
        "--top_file",
        help="this file provides topology information NOTE: when used with dcd must have .hoomdxml",
    )

    ##### Arguements for which frames to diffract. ######
    parser.add_argument(
        "-l",
        "--last_frame",
        help="generates diffraction data for the last frame of the traj file",
        action="store_true",
    )
    parser.add_argument(
        "-a",
        "--autocorr",
        help="this flag will use a autocorraltion function to generate diffraction data for each independent slice",
        action="store_true",
    )
    parser.add_argument(
        "-all", "--allf", help="diffract every attom frame", action="store_true"
    )

    ##### Arguements for which elements to include/exclude. ######
    parser.add_argument(
        "-ex",
        "--exclude_atom_name",
        help="exclude atoms listend with the -names/-name_file flag",
        action="store_true",
    )
    parser.add_argument(
        "-in",
        "--include_atom_name",
        help="include atoms listed with the -names/-names_file flag",
        action="store_true",
    )

    parser.add_argument(
        "-n",
        "--names",
        help="names of atoms to use in creating diffraction data, input names as a comma seperated list ie CHA,N,S or N",
    )
    parser.add_argument(
        "-nf",
        "--name_file",
        help="names of atoms to use in creating diffraction data, strcture file to have each atom name on a new line",
    )
    parser.add_argument(
        "-id",
        "--ids",
        help="ids of atoms to use in creating diffraction data, input ids as a comma seperated list ie 20,1,231 or 3",
    )
    parser.add_argument(
        "-idf",
        "--ids_file",
        help="ids of atoms to use in creating diffraction data, strcture file to have each atom ids on a new line",
    )

    ##### Arguements for naming the output directory. ######
    parser.add_argument(
        "-dir", "--dirname", help="Custom directory name. (Does not support spaces.)"
    )

    args = parser.parse_args()

    # This helps us reduce the number of cases

    atom_names = None
    atom_ids = None

    ##### For choosing atom types/ids ######

    if args.names is not None:
        args.names = args.names.split(",")  # formats args.names into a list
        atom_names = args.names

    if args.name_file is not None:
        atom_names = loadtxt(args.name_file, dtype=str)

    if args.ids is not None:
        args.ids = args.ids.split(",")  # formats args.names into a list
        atom_ids = args.ids

    if args.ids_file is not None:
        atom_ids = loadtxt(args.ids_file, dtype=int)

    ##### Calls to create diffraction data ######
    t = loadData(
        filetype=args.filetype, top_file=args.top_file, traj_file=args.traj_file
    )
    frames = generate_frames_array(
        last_frame=args.last_frame, autocorr=args.autocorr, allf=args.allf
    )
    iterate_through_desired_frames(
        t, frames, dirname=args.dirname, filetype=args.filetype
    )
