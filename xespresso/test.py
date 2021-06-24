import xml.etree.ElementTree as ET
from xmlschema import XMLSchema


xmlfile = '/home/xing/ase/xespresso/restart/scf/fe+u/fe+u.save/data-file-schema.xml'

xsd = XMLSchema('/home/xing/ase/xespresso/restart/scf/fe+u/fe+u.save/data-file-schema.xml')
# Validate XML document against the schema
# Returned dictionary has a structure where, if tag ['key'] is "simple", xml_dictionary['key'] returns its content.
# Otherwise, the following keys are available:
#
#  xml_dictionary['key']['$'] returns its content
#  xml_dictionary['key']['@attr'] returns its attribute 'attr'
#  xml_dictionary['key']['nested_key'] goes one level deeper.
xml_dictionary, errors = xsd.to_dict(xmlfile, validation='lax')
xml_version = xml_dictionary['general_info']['xml_format']['@VERSION']
inputs = xml_dictionary.get('input', {})
outputs = xml_dictionary['output']

print(inputs)
