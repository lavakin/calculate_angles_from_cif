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
from dint_golden_parameters import *
import json
import rmsd
import numpy as np
from array import array
from classes_for_angle_calculation import *


def calculate_euclidean_distance(dinc_tors: list, dinc_NCCN: list, step: list):
    dinc = [int(x) for x in dinc_tors] + [int(x) * 32 for x in dinc_NCCN]
    return round(np.linalg.norm(np.array(dinc) - np.array(step)),2)


ntcs = list(DINT_CLASS_INT.keys())[1:]


def GetAllAltpos(entity):
    variants = set()
    for sequence in entity.sequences:
        variants = variants.union(set(sequence.SequenceVariants.keys()))
    return variants


def GetNamesWithIndexes(entity):
    indexes = range(1, len(entity.steps) + 1)
    values = [x.name for x in entity.steps]
    return dict(zip(indexes, values))


def MakeJSON(models: [Model], meta_data, name):
    main_dict = {}
    with open("cif_data.json", "w") as json_file:
        names_with_indexes, euclid_distance = get_stepnrs_altpos_euclid(models)
        main_dict["stepnrs"] = names_with_indexes
        main_dict["eucl_dict"] = euclid_distance
        rmsd, C503, steps_amount = IterateTroughSteps(models)
        main_dict["pdbid"]= sys.argv[1][:-4]
        main_dict["rmsd_dict"] = rmsd
        main_dict["C503_dict"] = C503
        main_dict["prevnextids"] = get_prev_next(models)
        main_dict.update(steps_amount)
        main_dict["num_hetero_atoms"] = get_hetat(name)[0]
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
            for i in range(1,len(dic) + 1):
                for name in dic[i]:
                    prv_n_steps[name] = {}
                    if i != 1:
                        if len(dic[i]) == 1:
                            p = dic[i-1][0]
                        else:
                            p = next(x for x in dic[i-1] if parse_id_altpos(x)[1] == parse_id_altpos(name)[0])
                    else:
                        p = ""
                    prv_n_steps[name].update({'p':p})
                                             
                    if i != len(dic):
                        if len(dic[i]) == 1:
                            n = dic[i+1][0]
                        else:
                            n = next(x for x in dic[i+1] if parse_id_altpos(x)[0] == parse_id_altpos(name)[1])
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
    steps_amount = 0
    for model in models:
        for entity in model.entities:
            steps_amount += len(entity.steps)
            for step in entity.steps:
                seq_var1, seq_var2 = seq_id.find_seq_variants_by_step_id(models, step.name)
                rmsd[step.name], C503[step.name] = get_rmsd(step,seq_var1,seq_var2)
    return rmsd, C503, {"number_of_steps":steps_amount} 


def get_stepnrs_altpos_euclid(models):
    dict_of_euclid_for_model = {}
    names_with_indexes = {}
    for model in models:
        for entity in model.entities:
            names_with_indexes.update(GetNamesWithIndexes(entity))
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
