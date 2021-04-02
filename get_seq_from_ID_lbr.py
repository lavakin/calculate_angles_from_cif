#from calculate_angles_main import Name
import sys
from classes_for_angle_calculation import *

def find_seq_variants_by_step_id(models_array, name=None):
    model_number, name1, name2 = parse_step_ID(name)
    # first model, with matching number
    model = next(model for model in models_array if int(model.number) == model_number)
    seq_var1 = get_matched_seq_variant(model, name1)
    seq_var2 = get_matched_seq_variant(model, name2)
    return seq_var1, seq_var2


def get_matched_seq_variant(model, name):
    for entity in model.entities:
        for sequence in entity.sequences:
            # returns first seqence variant with the same name (equality redefined in the main program)
            sequence_var = \
                next((seq_var for seq_var in sequence.SequenceVariants.values() if seq_var.name == name), "none")
            if sequence_var != "none":
                return sequence_var
    return


def parse_step_ID(name=None):
    """
    parse the ID of the Step
    :return: model number and the name of both steps
    """
    if name==None:
        name=sys.argv[2]
    splitted = name.split("_")
    splitted_name = splitted[0].split("-")
    model_number = 1 if len(splitted_name) == 1 else int(splitted_name[1][1:])
    auth_asym_id = splitted[1]
    name1 = get_name_of_sequence(auth_asym_id, splitted[2], splitted[3])
    name2 = get_name_of_sequence(auth_asym_id, splitted[4], splitted[5])
    return model_number, name1, name2


def get_name_of_sequence(auth_asym_id, comp, seq_number):
    splitted_comp = comp.split(".")
    label_alt_id = "." if len(splitted_comp) == 1 else splitted_comp[1][:]
    splitted_seq_number = seq_number.split(".")
    inscode = "?" if len(splitted_seq_number) == 1 else splitted_seq_number[1][:]
    return Name(auth_asym_id, splitted_comp[0], label_alt_id, splitted_seq_number[0], inscode, "")



#_models_array = main.split_into_entities_and_calculate_parameters()
#_seq_var1, _seq_var2 = find_seq_variants_by_step_id()
