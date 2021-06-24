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
    """
    """
    for section, parameters in document.items():
        if key in parameters:
            parameter_type = document[section][key][0]
            # print(type(value))
            if parameter_type == "CHARACTER":
                   assert isinstance(value, str), '\n  Parameter: %s, should be a string!'%key
            elif parameter_type == 'REAL':
                   assert isinstance(value, float) or isinstance(value, int), '\n  Parameter: %s, should be a float!'%key
            elif parameter_type == 'INTEGER':
                   assert isinstance(value, int), '\n  Parameter: %s, should be a integer!'%key
            elif parameter_type == 'LOGICAL':
                   assert isinstance(value, bool), '\n  Parameter: %s, should be a bool!'%key
            option_list = document[section][key][1]
            if len(option_list) > 0:
                   assert value in option_list, '\n  Parameter: %s = %s, should be in %s!'%(key, value, str(option_list))
            return
    assert '\n  Parameter: %s, is not a PWSCF parameter!'%key
    
def modify_text(text, type = 'CHARACTER'):
    if text == None:
        text = ''
    else:
        text = str(text)
    text = text.replace(" ", "")
    text = text.replace("'", "").replace(":", "")
    text = text.strip()
    if type == 'LOGICAL' and text.upper() in '.TRUE.':
        text = True
    elif type == 'LOGICAL' and text.upper() in '.FALSE.':
        text = False
    if type.upper() == 'CHARACTER':
        if not text:
            text = ''
    if type.upper() == 'REAL':
        text = text.replace("D", "E")
        try:
            text = float(text)
        except:
            print('Failed: text %s'%text)
    if type.upper() == 'INTEGER':
        try:
            text = int(text)
        except:
            pass
    return text