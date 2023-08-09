def class_tree(cls):

    return { cls.__name__: [class_tree(sub_class) for sub_class in cls.__subclasses__()] }


def print_tree(tree, indent=4, current_ind=0):

    for k, v in tree.items():

        if current_ind:

            before_dashes = current_ind - indent
            print(' ' * before_dashes + '└' + '-'*(indent-1) + k)

        else:

            print(k)

        for sub_tree in v:

            print_tree(sub_tree, indent=indent, current_ind=current_ind + indent)


if __name__ == "__main__":

    from main_code.well_model.geometry_based_well_models.simple_pressure_losses_model import *
    from main_code.well_model.simplified_well.simplified_well import SimplifiedWell
    from main_code.well_model.geometry_based_well_models.REELWEEL_model import *
    from main_code.well_model.simplified_well.heating_sections import *

    print("------------------------- WELL MODELS -------------------------")
    ct = class_tree(SimplifiedWell)
    print_tree(ct, indent=4)
    print()
    print("---------------------- HEATING SECTIONS ----------------------")
    ct = class_tree(AbstractHeatingSection)
    print_tree(ct, indent=4)