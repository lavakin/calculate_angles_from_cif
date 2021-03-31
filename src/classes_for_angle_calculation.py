#!/usr/bin/python3
import sys
import math
from math import degrees
import gemmi
from gemmi import cif
import numpy as np
import csv
import locale
from copy import *
from locale import atof
import numpy

class Coords:
    def __init__(self, X, Y, Z):
        self.x = X
        self.y = Y
        self.z = Z



class Atom:
    def __init__(self, row):
        self.atom_name = remove_quotation_marks(row[3])
        self.coords = Coords(row[4], row[5], row[6])


class Model:
    def __init__(self, number):
        self.entities: [Entity] = []
        self.number = number


class Entity:
    def __init__(self):
        self.sequences: [Sequence] = []
        self.steps: [Step] = []


class Step:
    def __init__(self, model_number, is_the_only_one, seq1=None, seq2=None):
        self.name = get_step_name(seq1, seq2, model_number, is_the_only_one)
        self.angles = Angles(seq1, seq2)


class Angles:
    def __init__(self, seq1=None, seq2=None):
        """
        function computes all the angles for one step if there are enough data provided
        :param seq1: dict structure of the first sequence
        :param seq2: dict structure of the second sequence
        """
        if seq1 is None:
            self.has_angles = False
        else:
            self.has_angles = True
            self.alpha = calculate_dihedral_angles(seq1, seq2, ["O3'"], ["P", "O5'", "C5'"])
            self.beta = calculate_dihedral_angles(seq1, seq2, [], ["P", "O5'", "C5'", "C4'"])
            self.delta = calculate_dihedral_angles(seq1, seq2, ["C5'", "C4'", "C3'", "O3'"], [])
            self.epsilon = calculate_dihedral_angles(seq1, seq2, ["C4'", "C3'", "O3'"], ["P"])
            self.zeta = calculate_dihedral_angles(seq1, seq2, ["C3'", "O3'"], ["P", "O5'"])
            self.gamma = calculate_dihedral_angles(seq1, seq2, [], ["O5'", "C5'", "C4'", "C3'"])
            self.delta2 = calculate_dihedral_angles(seq1, seq2, [], ["C5'", "C4'", "C3'", "O3'"])
            self.chi = calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "N9", "C4"], []) \
                if seq1.base_type == "purine" \
                else calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "N1", "C2"], [])
            self.chi2 = calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "N9", "C4"]) \
                if seq2.base_type == "purine" \
                else calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "N1", "C2"])
            self.NCCN_tors = calculate_NCCN_torsion(seq1, seq2)
            self.taus = Taus(seq1, seq2)
            self.NNdist = get_NNdist(seq1, seq2)
            self.CCdist = get_distance(seq1.atoms["C1'"].coords, seq2.atoms["C1'"].coords)


class Taus:
    def __init__(self, seq1, seq2):
        self.tau0_1 = calculate_dihedral_angles(seq1, seq2, ["C4'", "O4'", "C1'", "C2'"], [])
        self.tau1_1 = calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "C2'", "C3'"], [])
        self.tau2_1 = calculate_dihedral_angles(seq1, seq2, ["C1'", "C2'", "C3'", "C4'"], [])
        self.tau3_1 = calculate_dihedral_angles(seq1, seq2, ["C2'", "C3'", "C4'", "O4'"], [])
        self.tau4_1 = calculate_dihedral_angles(seq1, seq2, ["C3'", "C4'", "O4'", "C1'"], [])
        self.tau0_2 = calculate_dihedral_angles(seq1, seq2, [], ["C4'", "O4'", "C1'", "C2'"])
        self.tau1_2 = calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "C2'", "C3'"])
        self.tau2_2 = calculate_dihedral_angles(seq1, seq2, [], ["C1'", "C2'", "C3'", "C4'"])
        self.tau3_2 = calculate_dihedral_angles(seq1, seq2, [], ["C2'", "C3'", "C4'", "O4'"])
        self.tau4_2 = calculate_dihedral_angles(seq1, seq2, [], ["C3'", "C4'", "O4'", "C1'"])
        self.tau1, self.t1_max = get_tau(self.tau0_1, self.tau1_1, self.tau2_1, self.tau3_1, self.tau4_1)
        self.tau1_type = get_pseudorotation_name(self.tau1)
        self.tau2, self.t2_max = get_tau(self.tau0_2, self.tau1_2, self.tau2_2, self.tau3_2, self.tau4_2)
        self.tau2_type = get_pseudorotation_name(self.tau2)


class Sequence:
    def __init__(self, asym_id, comp_id, seq_number, inscode, is_hetatom):
        self.SequenceVariants: {str: SequenceVariant} = {
            ".": SequenceVariant(asym_id, comp_id, ".", seq_number, inscode, is_hetatom)}


class Name:
    def __init__(self, asym_id, comp_id, alt_id, seq_number, inscode):
        self.seq_number = seq_number
        self.auth_asym_id = asym_id
        self.auth_comp_id = comp_id
        self.label_alt_id = alt_id
        self.inscode = inscode

    # redefining "==" to equality of properties
    def __eq__(self, other):
        if other is None:
            return False
        return (int(self.seq_number), self.auth_asym_id, self.auth_comp_id, self.label_alt_id, self.inscode) == \
               (int(other.seq_number), other.auth_asym_id, other.auth_comp_id, other.label_alt_id, other.inscode)


class SequenceVariant:
    def __init__(self, asym_id, comp_id, alt_id, seq_number, inscode, is_hetatom):
        self.atoms: {str: Atom} = {}
        self.is_valid = True
        self.base_type = None
        self.name = Name(asym_id, comp_id, alt_id, seq_number, inscode)
        self.is_hetatom = is_hetatom


class WrongNumberOfArgumentsException(Exception):
    def __init__(self, amount, message="Not the right amount of arguments"):
        self.amount = amount
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.amount} -> {self.message}'


def get_NNdist(seq1, seq2):
    """
    get the distance between two N on bases
    :param seq1:
    :param seq2:
    :return: euclidean distance
    """
    if seq1.base_type == "purine" and seq2.base_type == "purine":
        return get_distance(seq1.atoms["N9"].coords, seq2.atoms["N9"].coords)
    if seq1.base_type == "pyrimidine" and seq2.base_type == "purine":
        return get_distance(seq1.atoms["N1"].coords, seq2.atoms["N9"].coords)
    if seq1.base_type == "purine" and seq2.base_type == "pyrimidine":
        return get_distance(seq1.atoms["N9"].coords, seq2.atoms["N1"].coords)
    if seq1.base_type == "pyrimidine" and seq2.base_type == "pyrimidine":
        return get_distance(seq1.atoms["N1"].coords, seq2.atoms["N1"].coords)


def get_distance(atom1: Coords, atom2: Coords):
    """

    :param atom1: coords of the first atom
    :param atom2: coords of the second one
    :return: euclidean distance between two atoms
    """
    a = np.array((float(atom1.x), float(atom1.y), float(atom1.z)))
    b = np.array((float(atom2.x), float(atom2.y), float(atom2.z)))
    return round(np.linalg.norm(a - b), 2)


def get_tau(tau0, tau1, tau2, tau3, tau4):
    tau2 = 0.1 if tau2 == 0 else tau2
    tau0, tau1, tau2, tau3, tau4 = normalize_taus([tau0, tau1, tau2, tau3, tau4])
    tan_p = ((tau4 + tau1) - (tau3 + tau0)) / (
            2 * tau2 * (math.sin(36.0 * math.pi / 180.0) + math.sin(72.0 * math.pi / 180.0)))
    p = math.atan(tan_p) * 180.0 / math.pi
    if tau2 < 0:
        p = p + 180.0
    elif tan_p < 0:
        p = p + 360.0
    tmax = abs(tau2 / math.cos(p * math.pi / 180.0))
    return "%.1f" % p, "%.1f" % tmax


def normalize_taus(taus):
    taus = [tau - 360.0 if tau > 180.0 else tau for tau in taus]
    return taus


def get_pseudorotation_name(P):
    """
    For the given pseudorotation angle return its name. Based on Altona, Sundaralingam, JACS 94:23, 1972, 8205-8212
    :param P: pseudorotation
    :return: type of pseudorotation
    """
    P = atof(P)
    if ((P >= 0.0) and (P < 36.0)) or (P == 360.0):
        return "C3end"
    if (P >= 36.0) and (P < 72.0):
        return "C4exo"
    if (P >= 72.0) and (P < 108.0):
        return "O4end"
    if (P >= 108.0) and (P < 144.0):
        return "C1exo"
    if (P >= 144.0) and (P < 180.0):
        return "C2end"
    if (P >= 180.0) and (P < 216.0):
        return "C3exo"
    if (P >= 216.0) and (P < 252.0):
        return "C4end"
    if (P >= 252.0) and (P < 288.0):
        return "O4exo"
    if (P >= 288.0) and (P < 324.0):
        return "C1end"
    if (P >= 324.0) and (P < 360.0):
        return "C2exo"
    return ""


def get_step_name(seq1: SequenceVariant, seq2: SequenceVariant, model_number, is_the_only_one):
    seq1_alt = "." + seq1.name.label_alt_id if seq1.name.label_alt_id != "." else ""
    seq2_alt = "." + seq2.name.label_alt_id if seq2.name.label_alt_id != "." else ""
    model_num = "-m" + str(model_number) if not is_the_only_one else ""
    inscode1 = "." + seq1.name.inscode if seq1.name.inscode != "?" else ""
    inscode2 = "." + seq2.name.inscode if seq2.name.inscode != "?" else ""
    return (
            get_cif_name() + model_num + "_" + seq1.name.auth_asym_id + "_" + seq1.name.auth_comp_id + seq1_alt + "_" +
            seq1.name.seq_number + inscode1 + "_" + seq2.name.auth_comp_id + seq2_alt + "_" + seq2.name.seq_number
            + inscode2)


def calculate_NCCN_torsion(seq1, seq2):
    """
    calculates the NCCN torsion angle, which depends on purine or pyrimidine structure of sequences
    :param seq1: dict structure of the first sequence
    :param seq2: dict structure of the second sequence
    :return: a dihedral angle
    """
    if seq1.base_type == "purine" and seq2.base_type == "purine":
        return calculate_dihedral_angles(seq1, seq2, ["N9", "C1'"], ["C1'", "N9"])
    elif seq1.base_type == "pyrimidine" and seq2.base_type == "purine":
        return calculate_dihedral_angles(seq1, seq2, ["N1", "C1'"], ["C1'", "N9"])
    elif seq1.base_type == "purine" and seq2.base_type == "pyrimidine":
        return calculate_dihedral_angles(seq1, seq2, ["N9", "C1'"], ["C1'", "N1"])
    elif seq1.base_type == "pyrimidine" and seq2.base_type == "pyrimidine":
        return calculate_dihedral_angles(seq1, seq2, ["N1", "C1'"], ["C1'", "N1"])


def calculate_dihedral_angles(seq1: SequenceVariant, seq2: SequenceVariant, atom_names_seq1: [str],
                              atom_names_seq2: [str]):
    """
    calculates a dihedral angles of given points from two following sequences.
    The arrays may be empty, but it should be 4 points all together
    :param seq1: dict structure of the first sequence
    :param seq2: dict structure of the first sequence
    :param atom_names_seq1: the atoms that are needed from the first sequence
    :param atom_names_seq2: the atoms that are needed from the second sequence
    :return: an dihedral angle
    """

    points = [make_gemmi_position_format_from_coords(seq1.atoms[atom_name].coords) for atom_name in atom_names_seq1] + \
             [make_gemmi_position_format_from_coords(seq2.atoms[atom_name].coords) for atom_name in atom_names_seq2]
    if len(points) == 4:
        result = degrees(gemmi.calculate_dihedral(points[0], points[1], points[2], points[3]))
        return round(result, 1) if result >= 0 else round(360 + result, 1)  # angles should be > 0
    else:
        raise WrongNumberOfArgumentsException(len(points))


def remove_quotation_marks(atom):
    """
    solves the problem that cif stores the atoms with ' (ex. O4') with quotation marks, which needs to be removed
    :param atom: the name of the atom
    :return: name without quotation marks if there were any
    """
    return (atom[1:-1]) if '"' in atom else atom


def calculate_angles(entity: Entity, model_number, is_just_one_model):
    """
    Calculates angles if all needed atoms are present in the sequence, for all the sequence variants - for the same alt
    positions and for an alt position with canonical sequence.
    :param entity: one Entity
    :param model_number:
    :param is_just_one_model: for names sake (step ID)
    :return: an array of steps for the entity.
    """
    for i in range(len(entity.sequences) - 1):
        for sequence_var, sequence_var_object in entity.sequences[i].SequenceVariants.items():  # key and value
            if sequence_var_object.is_valid:
                if sequence_var == ".":  # if it is a canonical sequence
                    # calculates angles for all valid seq. variants of + 1 sequence.
                    for sequence_var2 in entity.sequences[i + 1].SequenceVariants.values():
                        if sequence_var2.is_valid:
                            entity.steps.append(
                                Step(model_number, is_just_one_model,
                                     sequence_var_object, sequence_var2))
                else:
                    # calculates angles for the matching alt position (if there is any valid one)
                    if sequence_var in entity.sequences[i + 1].SequenceVariants and \
                            entity.sequences[i + 1].SequenceVariants[sequence_var].is_valid:
                        entity.steps.append(Step(model_number, is_just_one_model, sequence_var_object,
                                                 entity.sequences[i + 1].SequenceVariants[sequence_var]))
                    # calculates angles for the canonical basis (if there is any valid one)
                    if "." in entity.sequences[i + 1].SequenceVariants and \
                            entity.sequences[i + 1].SequenceVariants["."].is_valid:
                        entity.steps.append(Step(model_number, is_just_one_model, sequence_var_object,
                                                 entity.sequences[i + 1].SequenceVariants["."]))


previous_comp_id = ""



def get_cif_name():
    # the cif can have ".cif" or ".gz" ending
    if sys.argv[1][-2:] != "gz":
        return sys.argv[1][:-4]
    else:
        return sys.argv[1][:-3]



def make_gemmi_position_format_from_coords(coords):
    return gemmi.Position(float(coords.x), float(coords.y), float(coords.z))


