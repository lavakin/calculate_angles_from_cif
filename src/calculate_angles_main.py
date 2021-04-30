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
#import store_to_json_cif as stjc
import get_seq_from_ID_lbr as seq_id
#sys.path.append("/srv/www/dnas/lib/backassign/pdbx")
from dint_golden_parameters import *
import json
import rmsd
import numpy
from store_to_json import *
from classes_for_angle_calculation import *
import classes_for_angle_calculation

purines = ["A","G","DA","DG","DDG","EDA","GNE","N2G","N5I","2DA","7DA","PRN","AD2","A3P","A5L","FMG","MA7","MG1","O2G","PPW","1AP","2FI","2PR","6MA","6OG","7GU","8OG","TGP","GFL","A2M","OMG","GTP","GDP","2MG","G7M","IGU","6NW","ATP"]
pyrimidines = ["T","C","U","DT","2DT","5NC","DC","DU","BRU","CBR","C38","DOC","ME6","OMC","UMP","Z","5CM","5IU","5PY","PST","SPT","TPC","TSP","UPS","US1","4PC","5HU","5FC","UFT","CFL","TAF","5HC","CCC","IMC","5BU","6OO","F2T","XFC"]


def parse_to_entities(table):
    """
    takes the whole structure and parses it to entities
    :param table: gemmi_table with all necessary information
    :return: array of entities
    """
    global previous_comp_id
    array_of_models = []
    numpy_table = np.array(table)
    m, n, k = "", "", ""
    for row in numpy_table:
        if row[10] != k:  # models
            array_of_models.append(Model(row[10]))
            k = row[10]
        if row[1] != n or (row[2] != m and row[2] == "1"):  # entities
            array_of_models[-1].entities.append(Entity())
            n = row[1]
            # append the first sequence of an entity
            array_of_models[-1].entities[-1].sequences.append(
                Sequence(row[8], row[9], row[11], row[12], row[0] == "HETATM", row[14]))
            m = row[2]
        if row[2] != m:  # sequences
            array_of_models[-1].entities[-1].sequences.append(
                Sequence(row[8], row[9], row[11], row[12], row[0] == "HETATM", row[14]))
            m = row[2]
        if row[7] not in array_of_models[-1].entities[-1].sequences[-1].SequenceVariants:  # new sequence variant
            array_of_models[-1].entities[-1].sequences[-1].SequenceVariants[row[7]] = \
                SequenceVariant(row[8], row[9], row[7], row[11], row[12], row[0] == "HETATM",row[14])
        # adding an atom into a sequence variant
        array_of_models[-1].entities[-1].sequences[-1].SequenceVariants[row[7]].atoms.update(
            {remove_quotation_marks(row[3]): Atom(row)})
    handle_altpos(array_of_models)
    return array_of_models


def handle_altpos(array_of_models):
    """
    removes void "." altpositions, adds the atoms from "." to the alt positions
    :param array_of_models: the whole structure
    :return:
    """
    for model in array_of_models:
        for entity in model.entities:
            for sequence in entity.sequences:
                if len(sequence.SequenceVariants) > 1:
                    without_altpos = copy(sequence.SequenceVariants["."])
                    del sequence.SequenceVariants["."]
                    for variant in sequence.SequenceVariants.values():
                        variant.atoms.update(without_altpos.atoms)


def operate_cif_file():
    doc = cif.read(sys.argv[1])
    block = doc[0]
    return block.find(
        ["_atom_site.group_PDB", "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.label_atom_id",
         "_atom_site.Cartn_x",
         "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.label_alt_id", "_atom_site.auth_asym_id",
         "_atom_site.auth_comp_id", "_atom_site.pdbx_PDB_model_num", "_atom_site.auth_seq_id",
         "_atom_site.pdbx_PDB_ins_code", "_atom_site.type_symbol","_atom_site.label_asym_id"])
    

def check_if_sequence_is_valid_and_add_basic_information(array_of_models):
    """
    the sequence is valid only if it contains all of the atoms given in variable "must_have". The function goes through
    all sequences in all entities and checks if sequences are valid. In addition to that it classifies the base as
    a purine or a pyrimidine.
    :param array_of_models: all entities in given structure
    :return: is_valid and base_type property of a sequence
    """
    must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "O4'", "C2'"]
    for model in array_of_models:
        for entity in model.entities:
            # the first sequence in entity does not have "P", that is why it is handled separately.
            for sequence_var in entity.sequences[0].SequenceVariants.values():
                all_atoms = sequence_var.atoms.keys()
                for atom in must_have:
                    if atom not in all_atoms and atom != "P":
                        sequence_var.is_valid = False
                if "N9" not in all_atoms and "N1" not in all_atoms:  # purine/pyrimidine validity
                    sequence_var.is_valid = False
                sequence_var.base_type = purine_or_pyrimidine(all_atoms)
                if sequence_var.is_hetatom:
                    # for hetatoms there is a special table with vald ones
                    if sequence_var.name.auth_comp_id not in pyrimidines and \
                            sequence_var.name.auth_comp_id not in purines:
                        sequence_var.is_valid = False

            for sequence in entity.sequences[1:]:
                for sequence_var in sequence.SequenceVariants.values():
                    all_atoms = [x for x in sequence_var.atoms.keys()]
                    for atom in must_have:
                        if atom not in all_atoms:
                            sequence_var.is_valid = False
                    if "N9" not in all_atoms and "N1" not in all_atoms:  # purine/pyrimidine validity
                        sequence_var.is_valid = False
                    sequence_var.base_type = purine_or_pyrimidine(all_atoms)
                    if sequence_var.is_hetatom:
                        if sequence_var.name.auth_comp_id not in pyrimidines and \
                                sequence_var.name.auth_comp_id not in purines:
                            sequence_var.is_valid = False


def purine_or_pyrimidine(all_atoms):
    if "N9" in all_atoms and "N1" in all_atoms:
        base_type = "purine"
    elif "N1" in all_atoms:
        base_type = "pyrimidine"
    else:
        base_type = "incomplete sequence"
    return base_type


def split_into_entities_and_calculate_parameters():
    cif_file = operate_cif_file()
    models_array = parse_to_entities(cif_file)
    meta_data = get_meta_data(cif_file)
    check_if_sequence_is_valid_and_add_basic_information(models_array)
    for model in models_array:
        for entity in model.entities:
            classes_for_angle_calculation.calculate_angles(entity, models_array.index(model) + 1, len(models_array) == 1)
    return models_array, meta_data


def write_to_csv():
    if not os.path.exists(os.getcwd() + '/my_csvs'):
        os.mkdir(os.getcwd() + "/my_csvs")
    with open(os.getcwd() + "/my_csvs/" + get_cif_name() + '.csv', 'w', newline='') as file:
        writer = csv.writer(file, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["step_ID", "d1", "e1", "z1", "a2", "b2", "g2", "d2", "ch1", "ch2", "NN", "CC", "mu", "P1",
                         "t1", "Pn1", "P2", "t2", "Pn2", "nu11", "nu12", "nu13", "nu14", "nu15", "nu21", "nu22", "nu23",
                         "nu24", "nu25"])
        for _model in _models_array:
            for _entity in _model.entities:
                for _atom in _entity.steps:
                    writer.writerow(
                        [_atom.name, _atom.angles.delta, _atom.angles.epsilon, _atom.angles.zeta, _atom.angles.alpha,
                         _atom.angles.beta,
                         _atom.angles.gamma, _atom.angles.delta2, _atom.angles.chi, _atom.angles.chi2,
                         _atom.angles.NNdist, _atom.angles.CCdist, _atom.angles.NCCN_tors, _atom.angles.taus.tau1,
                         _atom.angles.taus.t1_max, _atom.angles.taus.tau1_type, _atom.angles.taus.tau2,
                         _atom.angles.taus.t2_max, _atom.angles.taus.tau2_type, _atom.angles.taus.tau0_1,
                         _atom.angles.taus.tau1_1, _atom.angles.taus.tau2_1, _atom.angles.taus.tau3_1,
                         _atom.angles.taus.tau4_1, _atom.angles.taus.tau0_2,
                         _atom.angles.taus.tau1_2, _atom.angles.taus.tau2_2, _atom.angles.taus.tau3_2,
                         _atom.angles.taus.tau4_2])


def get_meta_data(cif_file):
    rows = np.transpose(cif_file)
    asym = {"chains":list(set([x[8] for x in np.transpose(rows) if len(x[9])==2]))} # to get rid of multiplicity
    alt_pos = {"alt_pos":list(set(rows[7]))}
    num_models = {"num_models" : len(set(rows[10]))}
    return [asym, alt_pos, num_models]


_models_array, meta_data = split_into_entities_and_calculate_parameters()
write_to_csv()
MakeJSON(_models_array, meta_data, sys.argv[1])
print("done")
