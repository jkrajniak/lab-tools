import collections
import logging
import networkx as nx
import pandas as pd
import re


logger = logging.getLogger()


class LammpsReader(object):
    """Very simple LAMMPS data file and input parser."""

    def __init__(self, verbose=True):
        self.verbose = verbose
        self.timestep = 0
        self.previous_section = None
        self.current_section = ""
        self._item_counters = {}
        self._type_counters = {}
        self._mass_type = {}
        self._section_line = None
        self.box = {}
        self.atoms = collections.defaultdict(dict)
        self.topology = {}
        self.distance_scale_factor = 0.1
        self.atom_charges = {}
        self.data_parsers = {
            "Atoms": self._read_atom,
            "Velocities": self._read_velocity,
            "Masses": self._read_mass,
            "Bonds": self._read_bond,
            "Angles": self._read_angle,
            "Dihedrals": self._read_dihedral,
            "Impropers": self._read_improper,
            "header": self._read_header,
            "coeffs": self._read_coeff,
        }
        self.force_field = collections.defaultdict(dict)
        self.init()
        # We shift box to the origin and we need to do the same with the atom position.
        self._box_translate = {}
        self.units = None

    def init(self):
        self.current_section = "header"
        self.previous_section = None

        # Data structures
        self._item_counters = {}
        self._type_counters = {}
        self._mass_type = {}
        self._section_line = None

        self.box = {}
        self.atoms = collections.defaultdict(dict)
        self.topology = {
            "bonds": collections.defaultdict(list),
            "angles": collections.defaultdict(list),
            "dihedrals": collections.defaultdict(list),
            "impropers": collections.defaultdict(list),
        }
        self.distance_scale_factor = 0.1

    def read_data(self, file_name, scale_factor=None, update=False):
        """Reads data file written with write_data command.

        Arsgs:
            file_name: The name of data file to read.
            scale_factor: The factor by which every distance quantity will be multiply.
        """
        if update:
            self.init()

        if scale_factor is not None:
            self.distance_scale_factor = scale_factor

        re_timestep = re.compile(".*timestep = ([0-9]+).*")

        with open(file_name, "r") as f:
            for line in f:
                line = line.strip().split("#")[0]
                if not line or line.startswith("#"):
                    continue
                if "timestep" in line:
                    timestep = re_timestep.match(line)
                    if timestep:
                        self.timestep = int(timestep.groups()[0])
                section_line = line.split("#")[0].strip()
                if section_line in self.data_parsers:
                    if self.verbose:
                        logger.info(("{}: Reading section {}".format(file_name, line)))
                    self.previous_section = self.current_section
                    self.current_section = section_line
                elif "Coeff" in section_line:
                    if self.verbose:
                        logger.info(("{}: Reading coefficient section {}".format(file_name, line)))
                    self.previous_section = self.current_section
                    self.current_section = "coeffs"
                    self._section_line = section_line
                elif self.current_section is not None:
                    self.data_parsers[self.current_section](line)
                else:
                    self.previous_section = self.current_section
                    self.current_section = None

    def print_info(self):
        """Prints information about the data."""
        logger.info("Timestep: %d", self.timestep)
        logger.info("Box: %s", self.box)
        logger.info("Atoms: %d", len(self.atoms))
        logger.info("Num atom types: %d", len(self._type_counters))
        logger.info("Topology:")
        logger.info(f"Num bond types: {len(self.topology['bonds'])}")
        logger.info(f"Num bonds: {sum([len(x) for x in self.topology['bonds'].values()])}")
        logger.info(f"Num angle types: {len(self.topology['angles'])}")
        logger.info(f"Num angles: {sum([len(x) for x in self.topology['angles'].values()])}")
        logger.info(f"Num dihedral types: {len(self.topology['dihedrals'])}")
        logger.info(f"Num dihedrals: {sum([len(x) for x in self.topology['dihedrals'].values()])}")
        logger.info(f"Num improper types: {len(self.topology['impropers'])}")
        logger.info(f"Num impropers: {sum([len(x) for x in self.topology['impropers'].values()])}")

    def read_dump(self, file_name, timestep, scale_factor=1.0, update=False):
        """Reads data file written with write_dump command.

        Args:
            file_name: The name of data file to read.
            timestep: The time step to read.
            scale_factor: The factor by which every distance is rescaled.
            update: If True then init() method is run first.
        """
        if update:
            self.init()

        if scale_factor is not None:
            self.distance_scale_factor = scale_factor

        current_item = None
        skip_frame = False
        with open(file_name, "r") as f:
            for line in f:
                if skip_frame and ("TIMESTEP" not in line and current_item != "TIMESTEP"):
                    continue
                if line.startswith("ITEM:"):
                    current_item = line.split(":")[1].strip()
                else:
                    if current_item == "TIMESTEP":
                        current_timestep = int(line)
                        if current_timestep != timestep:
                            skip_frame = True
                        else:
                            skip_frame = False
                    elif current_item == "NUMBER OF ATOMS":
                        int(line)
                    elif "ATOMS" in current_item:
                        atom_data = dict(list(zip(current_item.replace("ATOMS", "").split(), line.split())))
                        self.atoms[atom_data["id"]] = atom_data

    def read_input(self, file_name):
        """Reads LAMMPS input script. Only take cares on *_style and pair_coeff

        Args:
            file_name: Name of input file.
        """
        with open(file_name, "r") as f:
            for line in f:
                line = line.split("#")[0].strip()
                if not line or line.startswith("#"):  # skip empty lines
                    continue
                if "_style" in line:
                    sp_line = line.split()
                    if self.verbose:
                        logger.info(("Reads  {}".format(sp_line[0])))
                    self.force_field[sp_line[0]] = sp_line[1:]
                elif "bond_coeff" in line or "angle_coeff" in line or "dihedral_coeff" in line:
                    sp_line = line.split()
                    stype = sp_line[0].replace("_coeff", "")
                    btype = int(sp_line[1].strip())
                    self.force_field[stype][btype] = sp_line[2:]
                elif "pair_coeff" in line:
                    sp_line = line.split()
                    at_1, at_2 = sp_line[1:3]
                    if "*" not in sp_line[1]:
                        at_1 = int(at_1)
                    if "*" not in sp_line[2]:
                        at_2 = int(at_2)
                    self.force_field["pair_coeff"][tuple(sorted((at_1, at_2)))] = sp_line[3:]
                elif "units" in line:
                    self.units = line.split()[1]
                    if self.units == "real":
                        self.distance_scale_factor = 10**-1
                elif "read_data" in line:
                    sp_line = line.split()
                    data_file = sp_line[1].strip()
                    if self.verbose:
                        logger.info(("Reads data file: {}".format(data_file)))
                    self.read_data(data_file)

    def update_atoms(self, file_name):
        """Reads LAMMPS data file and updates only atom section.

        Args:
            file_name: Input data file with Atoms section.
        """
        with open(file_name, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                section_line = line.split("#")[0].strip()
                if section_line in self.data_parsers:
                    self.previous_section = self.current_section
                    self.current_section = section_line
                    if section_line == "Atoms":
                        if self.verbose:
                            logger.info(('{}: Found "Atoms" section, updating atoms'.format(file_name)))
                elif self.current_section is not None and self.current_section == "Atoms":
                    self._read_atom(line, update=True)
                else:
                    self.current_section = None

    def get_atoms_as_dataframe(self) -> pd.DataFrame:
        """Returns atoms as pandas DataFrame."""
        return pd.DataFrame(self.atoms)

    def get_graph(self, settings):
        """Creates nx.Graph object from coordinate and topology data.

        Args:
            settings: The settings object.

        Returns:
            nx.Graph object. Each of node has attributes:
                - name: The name of atom.
                - res_id: The id of molecule.
                - chain_name: The name of molecule.
                - position: The x, y, z position tuple.
                - degree: The degree of the node.
            Graph object has attribute `box` with the x, y, z values of a box.
            Each of edges has attribute `bond_type` with a number that corresponds to bond type
            in force field data.
        """
        type2chain_name = settings.type2chain
        name_seq = settings.name_seq
        output_graph = nx.Graph(box=(self.box["x"], self.box["y"], self.box["z"]))
        seq_idx = {k: 0 for k in name_seq}
        for at_id, lmp_at in list(self.atoms.items()):
            chain_name = type2chain_name[lmp_at["atom_type"]]
            at_seq = name_seq[chain_name]
            chain_len = len(at_seq)
            at_name = at_seq[seq_idx[chain_name] % chain_len]
            mol_idx = lmp_at["res_id"]
            position = lmp_at["position"]
            output_graph.add_node(at_id, name=at_name, res_id=mol_idx, position=position, chain_name=chain_name)
            seq_idx[chain_name] += 1

        # Adding edges
        for bond_id, bond_list in list(self.topology["bonds"].items()):
            for b1, b2 in bond_list:
                output_graph.add_edge(b1, b2, bond_type=bond_id)

        # Updates degree
        for n_id in output_graph.nodes():
            output_graph.node[n_id]["degree"] = output_graph.degree(n_id)

        return output_graph

    def get_simple_graph(self):
        output_graph = nx.Graph(box=(self.box["x"], self.box["y"], self.box["z"]))
        for at_id, lmp_at in list(self.atoms.items()):
            atom_type = lmp_at["atom_type"]
            mol_idx = lmp_at["res_id"]
            position = lmp_at["position"]
            output_graph.add_node(at_id, position=position, atom_type=atom_type, mol_idx=mol_idx)
        # Adding edges
        for bond_id, bond_list in list(self.topology["bonds"].items()):
            for b1, b2 in bond_list:
                output_graph.add_edge(b1, b2, bond_type=bond_id)

        # Updates degree
        for n_id in output_graph.nodes():
            output_graph.nodes[n_id]["degree"] = output_graph.degree(n_id)

        return output_graph

    def write(self, output: str):
        g = self.get_simple_graph()
        content = ["Polymer Network System", ""]

        num_bonds = len([p for v in list(self.topology["bonds"].values()) for p in v])
        num_angles = len([p for v in list(self.topology["angles"].values()) for p in v])
        num_atom_types = max(set([a["atom_type"] for a in list(self.atoms.values())]))
        num_bond_types = len(self.topology["bonds"])
        num_angle_types = len(self.topology["angles"])

        content.append(f"{len(g.nodes)} atoms")
        content.append(f"{num_bonds} bonds")
        content.append(f"{num_angles} angles")
        content.append("")
        content.append(f"{num_atom_types} atom types")
        content.append(f"{num_bond_types} bond types")
        content.append(f"{num_angle_types} angle types")
        content.append("")
        content.append("0.0000000 {} xlo xhi".format(self.box["x"] / self.distance_scale_factor))
        content.append("0.0000000 {} ylo yhi".format(self.box["y"] / self.distance_scale_factor))
        content.append("0.0000000 {} zlo zhi".format(self.box["z"] / self.distance_scale_factor))
        content.append("")
        content.append("Masses")
        content.append("")
        for at_type, mass in list(self._mass_type.items()):
            content.append(f"{at_type} {mass}")
        content.append("")
        content.append("Atoms")
        content.append("")
        atom_ids = sorted(self.atoms.keys())
        for at_id in atom_ids:
            at_data = self.atoms[at_id]
            content.append(
                "{:8d} {:8d} {:8d} {:10.6f} {:8.3f} {:8.3f} {:8.3f}".format(
                    at_id,
                    at_id,
                    at_data["atom_type"],
                    0.0,
                    at_data["position"][0] / self.distance_scale_factor,
                    at_data["position"][1] / self.distance_scale_factor,
                    at_data["position"][2] / self.distance_scale_factor,
                )
            )
        content.append("")
        content.append("Bonds")
        content.append("")
        bond_idx = 1
        for bond_id, bond_list in list(self.topology["bonds"].items()):
            for b1, b2 in bond_list:
                content.append("{:8d} {:8d} {:8d} {:8d}".format(bond_idx, bond_id, b1, b2))
                bond_idx += 1
        content.append("")
        content.append("Angles")
        content.append("")
        angle_idx = 1
        for angle_id, angle_list in list(self.topology["angles"].items()):
            for a1, a2, a3 in angle_list:
                content.append("{:8d} {:8d} {:8d} {:8d} {:8d}".format(angle_idx, angle_id, a1, a2, a3))
                angle_idx += 1
        content.append("")

        with open(output, "w") as f:
            f.write("\n".join(content))

    # Parsers section
    def _read_header(self, input_line):
        """Parses header of data file."""
        sp_line = input_line.split()
        if "types" in sp_line:
            self._type_counters[sp_line[1]] = int(sp_line[0])
        elif "xhi" in sp_line or "yhi" in sp_line or "zhi" in sp_line:
            lo, hi = list(map(float, sp_line[:2]))
            lo *= self.distance_scale_factor
            hi *= self.distance_scale_factor
            tag = sp_line[-1].replace("hi", "")
            self._box_translate[tag] = lo
            self.box[tag] = hi - lo
        elif (
            "atoms" in sp_line
            or "bonds" in sp_line
            or "angles" in sp_line
            or "dihedrals" in sp_line
            or "impropers" in sp_line
        ):
            self._item_counters[sp_line[1]] = int(sp_line[0])

    def _read_coeff(self, input_line):
        coeff_type, _ = self._section_line.split()
        coeff_type = coeff_type.lower()

        sp_line = input_line.split()
        ff_type = int(sp_line[0])
        self.force_field[coeff_type][ff_type] = sp_line[1:]

    def _read_atom(self, input_line, update=False):
        sp_line = input_line.split()
        # Set type
        sp_line[:3] = list(map(int, sp_line[:3]))
        sp_line[3:7] = list(map(float, sp_line[3:7]))
        sp_line_len = len(sp_line)
        if sp_line_len == 10:
            sp_line[7:10] = list(map(int, sp_line[7:10]))
            at_id, at_tag, at_type, q, x, y, z, nx, ny, nz = sp_line
        elif sp_line_len == 6:
            at_id, at_tag, at_type, x, y, z = sp_line
            q = 0
            nx, ny, nz = None, None, None
        else:
            at_id, at_tag, at_type, q, x, y, z = sp_line
            nx, ny, nz = None, None, None

        if at_id > self._item_counters["atoms"]:
            raise RuntimeError(
                ('Number of atoms in "header" section does not ' 'correspond to number of atoms in "Atoms" section.')
            )

        if at_type > self._type_counters["atom"]:
            raise RuntimeError(("Atom type {} not found.".format(at_type)))

        # if at_type in self.atom_charges:
        #    if q != self.atom_charges[at_type]:
        #        raise RuntimeError('Charge of atom type {} is different, {} != {}'.format(
        #            at_type, q, self.atom_charges[at_type]
        #        ))
        self.atom_charges[at_type] = q

        if update:  # Update
            if at_id not in self.atoms:
                raise RuntimeError("Cannot update atom with id {}. Not found.".format(at_id))
            update_dict = {
                "position": (
                    x * self.distance_scale_factor,  # - self._box_translate['x'],
                    y * self.distance_scale_factor,  # - self._box_translate['y'],
                    z * self.distance_scale_factor,
                ),  # - self._box_translate['z']),
                "atom_type": at_type,
                "res_id": at_tag,
                "charge": q,
            }
            if nx is not None:
                update_dict["image"] = (nx, ny, nz)
            self.atoms[at_id].update(update_dict)
        else:  # New entry
            if at_id in self.atoms:
                raise RuntimeError("Cannot overwrite atom with id {} if update=False".format(at_id))
            self.atoms[at_id] = {
                "atom_type": at_type,
                "res_id": at_tag,
                "position": (
                    x * self.distance_scale_factor,  # - self._box_translate['x'],
                    y * self.distance_scale_factor,  # - self._box_translate['y'],
                    z * self.distance_scale_factor,
                ),  # - self._box_translate['z']),
                "image": (nx, ny, nz),
                "charge": q,
                "vel": (0.0, 0.0, 0.0),
                "mass": self._mass_type.get(at_type, 0.0),
            }

    def _read_velocity(self, input_line):
        sp_line = input_line.split()
        sp_line[0] = int(sp_line[0])
        sp_line[1:] = list(map(float, sp_line[1:]))
        at_id, vx, vy, vz = sp_line
        self.atoms[at_id]["vel"] = (
            vx * self.distance_scale_factor,
            vy * self.distance_scale_factor,
            vz * self.distance_scale_factor,
        )

    def _read_bond(self, input_line):
        idd, bond_type, at_1, at_2 = list(map(int, input_line.split()))
        if idd > self._item_counters["bonds"]:
            raise RuntimeError("Number of bond is wrong.")

        if at_1 not in self.atoms or at_2 not in self.atoms:
            raise RuntimeError("{} or {} not found in list of atoms.".format(at_1, at_2))
        self.topology["bonds"][bond_type].append(tuple(sorted((at_1, at_2))))

    def _read_angle(self, input_line):
        idd, angle_type, at_1, at_2, at_3 = list(map(int, input_line.split()))
        if idd > self._item_counters["angles"]:
            raise RuntimeError("Number of angle is wrong.")

        if at_1 not in self.atoms or at_2 not in self.atoms or at_3 not in self.atoms:
            raise RuntimeError("{}, {} or {} not found in list of atoms.".format(at_1, at_2, at_3))
        self.topology["angles"][angle_type].append((at_1, at_2, at_3))

    def _read_dihedral(self, input_line):
        idd, dihedral_type, at_1, at_2, at_3, at_4 = list(map(int, input_line.split()))
        if idd > self._item_counters["dihedrals"]:
            raise RuntimeError("Number of dihedrals is wrong.")

        if at_1 not in self.atoms or at_2 not in self.atoms or at_3 not in self.atoms or at_4 not in self.atoms:
            raise RuntimeError("{}, {}, {} or {} not found in list of atoms.".format(at_1, at_2, at_3, at_4))
        self.topology["dihedrals"][dihedral_type].append((at_1, at_2, at_3, at_4))

    def _read_improper(self, input_line):
        idd, dihedral_type, at_1, at_2, at_3, at_4 = list(map(int, input_line.split()))
        if idd > self._item_counters["impropers"]:
            raise RuntimeError("Number of impropers is wrong.")

        if at_1 not in self.atoms or at_2 not in self.atoms or at_3 not in self.atoms or at_4 not in self.atoms:
            raise RuntimeError("{}, {}, {} or {} not found in list of atoms.".format(at_1, at_2, at_3, at_4))
        self.topology["impropers"][dihedral_type].append((at_1, at_2, at_3, at_4))

    def _read_mass(self, input_line):
        sp_line = input_line.split()
        at_id, mass = sp_line
        self._mass_type[int(at_id)] = float(mass)
