import numpy as np


def get_hash(file):
    import hashlib

    with open(file, "rb") as f:
        file_hash = hashlib.md5()
        chunk = f.read(8192)
        while chunk:
            file_hash.update(chunk)
            chunk = f.read(8192)
    hd = file_hash.hexdigest()
    return hd


def check_type(key, value, document):
    """ """
    for section, parameters in document.items():
        if key in parameters:
            parameter_type = document[section][key][0]
            # print(type(value))
            if parameter_type == "CHARACTER":
                assert isinstance(value, str), (
                    "\n  Parameter: %s, should be a string!" % key
                )
            elif parameter_type == "REAL":
                assert isinstance(value, float) or isinstance(value, int), (
                    "\n  Parameter: %s, should be a float!" % key
                )
            elif parameter_type == "INTEGER":
                assert isinstance(value, int), (
                    "\n  Parameter: %s, should be a integer!" % key
                )
            elif parameter_type == "LOGICAL":
                assert isinstance(value, bool), (
                    "\n  Parameter: %s, should be a bool!" % key
                )
            option_list = document[section][key][1]
            if len(option_list) > 0:
                assert (
                    value in option_list
                ), "\n  Parameter: %s = %s, should be in %s!" % (
                    key,
                    value,
                    str(option_list),
                )
            return
    assert "\n  Parameter: %s, is not a PWSCF parameter!" % key


def modify_text(text, type="CHARACTER"):
    if text == None:
        text = ""
    else:
        text = str(text)
    text = text.replace(" ", "")
    text = text.replace("'", "").replace(":", "")
    text = text.strip()
    if type == "LOGICAL" and text.upper() in ".TRUE.":
        text = True
    elif type == "LOGICAL" and text.upper() in ".FALSE.":
        text = False
    if type.upper() == "CHARACTER":
        if not text:
            text = ""
    if type.upper() == "REAL":
        text = text.replace("D", "E")
        try:
            text = float(text)
        except:
            print("Failed: text %s" % text)
    if type.upper() == "INTEGER":
        try:
            text = int(text)
        except:
            pass
    return text


def compare_value(v1, v2, tol=1e-5):
    """Compare two value."""
    if isinstance(v1, str):
        if v1.upper() != v2.upper():
            return False
    elif isinstance(v1, dict):
        changed_parameters, igonre_parameters = compare_dict(v1, v2)
        if changed_parameters:
            return False
        for key, value in v1.items():
            if not compare_value(value, v2[key]):
                return False
    elif isinstance(v1, bool):
        if v1 != v2:
            return False
    else:
        if abs(v1 - v2) > tol:
            return False
    return True


def compare_dict(dict1, dict2, ignore=[], default=None):
    """Compare two dict
    Args:
        dict1 (_type_): _description_
        dict2 (_type_): _description_
        ignore (list, optional): _description_. Defaults to [].
        default (_type_, optional): _description_. Defaults to None.
    Returns:
        _type_: _description_
    """
    igonre_parameters = []
    changed_parameters = []
    keys = set(list(dict1.keys()) + list(dict2.keys()))
    for key in keys:
        if key in ignore:
            igonre_parameters.append(key)
        elif key not in dict1:
            changed_parameters.append(key)
        elif key not in dict2:
            if default and compare_value(dict1[key], default[key]):
                continue
            changed_parameters.append(key)
        elif not compare_value(dict1[key], dict2[key]):
            changed_parameters.append(key)
    return changed_parameters, igonre_parameters


def compare_parameters(para1, para2, ignore=[]):
    """ """
    from xespresso.input_parameters import (
        qe_namespace,
        default_parameters,
        restart_ignore,
    )

    changed_parameters = []
    igonre_parameters = []
    if not para1:
        changed_parameters = ["all"]
        return changed_parameters, igonre_parameters
    default_parameters = default_parameters["PW"]
    # pseudopotentials
    key = "pseudopotentials"
    try:
        for species, value in para1[key].items():
            if value != para2[key][species]:
                changed_parameters.append(key)
                continue
    except:
        changed_parameters.append(key)
    # kpts
    key = "kpts"
    try:
        if para1[key] != para2[key]:
            changed_parameters.append(key)
    except:
        changed_parameters.append(key)
    # input_data
    for section, paras in para1["input_data"].items():
        if section == "INPUT_NTYP":
            changed_parameters1, igonre_parameters1 = compare_dict(
                para1["input_data"][section],
                para2["input_data"][section],
                restart_ignore["PW"],
            )
        else:
            changed_parameters1, igonre_parameters1 = compare_dict(
                para1["input_data"][section],
                para2["input_data"][section],
                restart_ignore["PW"],
                default=default_parameters[section],
            )
        changed_parameters.extend(changed_parameters1)
        igonre_parameters.extend(igonre_parameters1)
    return changed_parameters, igonre_parameters
