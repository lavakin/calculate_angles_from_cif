#!/usr/bin/python3
import sys
import math
from math import degrees
from gemmi import cif
import numpy as np
import csv
from copy import *
import get_seq_from_ID_lbr as seq_id
from dint_golden_parameters import *
import json
import rmsd
from array import array
from classes_for_angle_calculation import *
import os



def periodic_distance(angle, angle2):
    return min(abs(angle2 - angle), 360 - (abs(angle2 - angle)))

def calculate_euclidean_distance(dinc_tors: list, dinc_NCCN: list, step: list):
    p = [periodic_distance(x,y) for (x,y) in zip(dinc_tors, step)] + [(abs(x-y)*32) for (x,y) in zip(dinc_NCCN[:-1], step[-3:-1])] + [periodic_distance(dinc_NCCN[-1],step[-1])]
    p = math.sqrt(sum(x**2 for x in p))
    return round(p,2)
    #return round(np.linalg.norm(np.array(dinc) - np.array(step)),2)

"""
def euclid():
    for torsidx, torsion in enumerate(torsions[0:9]):
    euclid += (distance_two_angles(float(data_matrix[num_steps - 1][torsion]),
                                    DINT_BBTORS_MEAN[DINT_CLASS_INT[ntc][0]][torsidx])) ** 2.0
    for torsidx, torsion in enumerate(torsions[9:12]):
        if (torsidx < 2):
            euclid += (dist_multiplier * distance_two_angles(
                float(data_matrix[num_steps - 1][torsion]),
                DINT_NCCN_MEAN[DINT_CLASS_INT[ntc][0]][torsidx])) ** 2.0
        else:
            euclid += (distance_two_angles(float(data_matrix[num_steps - 1][torsion]),
                                            DINT_NCCN_MEAN[DINT_CLASS_INT[ntc][0]][torsidx])) ** 2.0
    eucl_dict[step_ID][ntc] = "%.1f" % (math.sqrt(euclid))
"""


ntcs = list(DINT_CLASS_INT.keys())[1:]

def remove_quotation_marks(atom):
    """
    solves the problem that cif stores the atoms with ' (ex. O4') with quotation marks, which needs to be removed
    :param atom: the name of the atom
    :return: name without quotation marks if there were any
    """
    return atom.replace('"',"").replace("'","")

def GetAllAltpos(entity):
    variants = set()
    for sequence in entity.sequences:
        variants = variants.union(set(sequence.SequenceVariants.keys()))
    return variants


def GetNamesWithIndexes(entity, start_index):
    indexes = range(start_index, len(entity.steps) + start_index)
    values = [x.name for x in entity.steps]
    return dict(zip(indexes, values))


def MakeJSON(models: [Model], meta_data, name):
    print(os.getcwd())
    if not os.path.exists(os.getcwd()+'/my_jsons'):
        os.mkdir(os.getcwd()+"/my_jsons")
    doc = cif.read(name)
    block = doc[0]
    main_dict = {}
    title = np.array(block.find(["_entity.pdbx_description"]))[0][0][1:-1]
    het_at = get_proteins_and_het()

    with open(os.getcwd() + "/my_jsons/"+ get_cif_name()+".json", "w") as json_file:
        names_with_indexes, euclid_distance = get_stepnrs_altpos_euclid(models)
        main_dict["stepnrs"] = names_with_indexes
        main_dict["eucl_dict"] = euclid_distance
        rmsd, C503, steps_amount, properties = IterateTroughSteps(models)
        main_dict["pdbid"]= get_cif_name()
        main_dict["rmsd_dict"] = rmsd
        main_dict["C503_dict"] = C503
        main_dict["struct_title"] = str(title)
        main_dict["steps"] = properties
        main_dict["prevnextids"] = get_prev_next(models)
        main_dict.update(steps_amount)
        print(meta_data)
        #main_dict["num_hetero_atoms"] = get_hetat(name)[0]
        main_dict.update(het_at)
        print(het_at)
        main_dict["number_of_steps"] = len(names_with_indexes)
        for dictt in meta_data:
            main_dict.update(dictt)
        json.dump(main_dict, json_file)


def get_hetat(name):
    doc = cif.read(name)
    block = doc[0]
    m = np.array(block.find(["_entity.type"]))
    het =  len([x for x in m if x in  ["non-polymer", "macrolide"]])
    water = len([x for x in m if x == "water"])
    return het, water

def get_prev_next(models):
    prv_n_steps = {}
    for model in models:
        for entity in model.entities:
            dic = {}
            for i in range(len(entity.steps)):
                num = parse_id_seq(entity.steps[i].name)
                try:
                    dic[num].append(entity.steps[i].name)
                except:
                    dic[num] = [entity.steps[i].name]
            keys = list(dic.keys())
            for i in keys:
                for name in dic[i]:
                    prv_n_steps[name] = {}
                    if keys.index(i) != 0:
                        if len(dic[i]) == 1:
                            p = dic[keys[keys.index(i)-1]][0]
                        else:
                            p = next(x for x in dic[keys[keys.index(i)-1]] if parse_id_altpos(x)[1] == parse_id_altpos(name)[0] or len(dic[keys[keys.index(i)-1]])==1)
                    else:
                        p = ""
                    prv_n_steps[name].update({'p':p})
                                             
                    if keys.index(i) < len(dic)-1:
                        if len(dic[i]) == 1:
                            n = dic[keys[keys.index(i)+1]][0]
                        else:
                            
                            n = next(x for x in dic[keys[keys.index(i)+1]] if parse_id_altpos(x)[0] == parse_id_altpos(name)[1] or len(dic[keys[keys.index(i)+1]])==1)
                        
                    else:
                        n = ""
                    prv_n_steps[name].update({'n':n})
    return prv_n_steps
                    
                    
                
def parse_id_altpos(id):
    splitted = id.split("_")
    s1 = '.' if len(splitted[2].split("."))==1 else splitted[2].split(".")[1]
    s2 = '.' if len(splitted[4].split("."))==1 else splitted[4].split(".")[1]
    return s1,s2


def parse_id_seq(id):
    splitted = id.split("_")
    return int(splitted[3].split(".")[0])


def IterateTroughSteps(models):
    rmsd = {}
    C503 = {}
    steps = {}
    steps_amount = 0
    for model in models:
        for entity in model.entities:
            steps_amount += len(entity.steps)
            for step in entity.steps:
                seq_var1, seq_var2 = seq_id.find_seq_variants_by_step_id(models, step.name)
                properties = get_properties(seq_var1,seq_var2)
                steps[step.name] = properties 
                rmsd[step.name], C503[step.name] = get_rmsd(step,seq_var1,seq_var2)
    return rmsd, C503, {"number_of_steps":steps_amount}, steps


def get_properties(seq1,seq2):
    properties = {'label_alt_id_1':seq1.name.label_alt_id,'label_comp_id_1':seq1.name.auth_comp_id, 'PDB_ins_code_1':seq1.name.inscode, 'label_seq_id_1': seq1.name.seq_number, 'label_asym_id_1': seq1.name.label_asym_id, 'auth_asym_id_1':seq1.name.auth_asym_id, 'label_alt_id_2': seq2.name.label_alt_id,
                  'label_comp_id_2': seq2.name.auth_comp_id, 'PDB_ins_code_2': seq2.name.inscode, 'label_seq_id_2': seq2.name.seq_number, 'label_asym_id_2': seq2.name.label_asym_id, 'auth_asym_id_2': seq2.name.auth_asym_id }
    return properties
    
    
def get_stepnrs_altpos_euclid(models):
    dict_of_euclid_for_model = {}
    names_with_indexes = {}
    for model in models:
        for entity in model.entities:
            names_with_indexes.update(GetNamesWithIndexes(entity,max(names_with_indexes.keys(),default=0)+1))
            for step in entity.steps:
                dict_of_euclid = {}
                for index in DINT_NCCN_MEAN.keys():
                    dict_of_euclid[DINT_CLASS_NTC_CANA[index][0]] = calculate_euclidean_distance(
                        DINT_BBTORS_MEAN[index], DINT_NCCN_MEAN[index],
                        [step.angles.delta, step.angles.epsilon, step.angles.zeta, step.angles.alpha,
                         step.angles.beta,
                         step.angles.gamma, step.angles.delta2, step.angles.chi, step.angles.chi2,
                         step.angles.NNdist, step.angles.CCdist, step.angles.NCCN_tors])
                        
                dict_of_euclid_for_model[step.name] = dict_of_euclid
    return  names_with_indexes, dict_of_euclid_for_model

def get_rmsd(step, seq_var1, seq_var2):
    atoms1, atoms2 = get_coords(step,seq_var1,seq_var2)
    atoms = [atoms1[x].coords for x in atoms1.keys()] + [atoms2[x].coords for x in atoms2.keys()]
    atoms_in_vectors = [get_vectors_from_coords(x) for x in atoms]
    return calculate_rmsd(atoms_in_vectors)

def calculate_rmsd(atoms):
    rmsd_dict = {}
    C_dict = {}
    for ntc in ntcs:
        S = numpy.array(DINT_CLASS_COORDS[ntc])
        R = numpy.array(atoms)
        SC = rmsd.centroid(S)
        S -= SC
        RC = rmsd.centroid(R)
        R -= RC
        S = rmsd.kabsch_rotate(S, R)
        rmsd_dict[ntc] = "%.3f" % rmsd.rmsd(S, R)
        S += RC
        C_dict[ntc] = [["{:0.3f}".format(x) for x in S[0].tolist()],
                                   ["{:0.3f}".format(x) for x in S[3].tolist()],
                                   ["{:0.3f}".format(x) for x in S[10].tolist()],
                                   ["{:0.3f}".format(x) for x in S[13].tolist()]]
        
    return rmsd_dict, C_dict

def get_vectors_from_coords(coords: Coords):
    return [float(coords.x), float(coords.y), float(coords.z)]

def get_coords(step, seq_var1, seq_var2):
    atoms1 = {x:seq_var1.atoms[x] for x in ["C5'", "C4'", "C3'", "O3'", "O4'", "C1'"]}
    atoms2 = {x:seq_var2.atoms[x] for x in [ "P", "O5'", "C5'", "C4'", "C3'", "O3'", "O4'", "C1'" ]}
    atoms1 = add_pur_pyr(seq_var1,atoms1)
    atoms2 = add_pur_pyr(seq_var2,atoms2)
    return atoms1,atoms2

def add_pur_pyr(seq_variant, atoms):
    if seq_variant.base_type == 'purine':
        atoms.update({x:seq_variant.atoms[x] for x in ["N9", "C4" ]})
    else: atoms.update({x:seq_variant.atoms[x] for x in [ "N1", "C2" ]})
    return atoms


def get_proteins_and_het():
    entities_type = {'num_nucleic_atoms': 0, 'num_protein_atoms': 0, 'num_hetero_atoms' : 0, 'num_water_atoms': 0}
    # cif reader
    doc = cif.read(sys.argv[1])
    block = doc[0]
    entities = np.array(block.find(
        ["_atom_site.label_entity_id"]))
    # convert all values to int 
    entities = [int(x) for x in entities]
    # get entities in for of [(entity, number of atoms)]
    entities =  [ (i,entities.count(i)) for i in set(entities) ]
    properties = np.transpose(np.array(block.find(["_entity_poly.entity_id", "_entity_poly.type"])))
    properties[1] = [remove_quotation_marks(x) for x in properties[1]]
    print(properties)
    properties2 = np.transpose(np.array(block.find(["_entity.id", "_entity.type"])))
    # make dict, with entity number as a key and property as a value
    properties = dict(zip([int(x) for x in properties[0]], [[x] for x in properties[1]]))
    properties2 = dict(zip([int(x) for x in properties2[0]], [[x] for x in properties2[1]]))

    # merge properties in one dict - value is a list of properties for given entity
    for index,proper in properties2.items():
        if index in properties.keys():
            properties[index] = properties[index] + proper
        else:
            properties[index] = proper

    #count number of atoms belonging to each group
    for entity in entities:
        # set "and" operator - true if 2 lists have at least one element in common
        if (set(properties[entity[0]]) & set(["polydeoxyribonucleotide", "polyribonucleotide", "polydeoxyribonucleotide/polyribonucleotide hybrid"])):
            entities_type['num_nucleic_atoms'] += entity[1]
        elif (set(properties[entity[0]]) & set(["polypeptide(D)", "polypeptide(L)"])):
            entities_type['num_protein_atoms'] += entity[1]
        elif set(properties[entity[0]]) & set(["polysaccharide(D)", "polysaccharide(L)", "cyclic-pseudo-peptide", "other", "peptide nucleic acid", "non-polymer", "macrolide"]):
            entities_type['num_hetero_atoms'] += entity[1]
        elif ("water" in properties[entity[0]]):
            entities_type['num_water_atoms'] += entity[1]
    print(entities_type)
    return entities_type
    

    
