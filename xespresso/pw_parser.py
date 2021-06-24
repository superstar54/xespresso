"""
Creating a json file from the input_pw.html file
"""

import requests
from bs4 import BeautifulSoup
from requests.api import options
import json
# from xespresso.utils import modify_text

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
            print('Failed:   %s    %s'%(text, type))
            text = 1e30
    if type.upper() == 'INTEGER':
        try:
            text = int(text)
        except:
            print('Failed:   %s    %s'%(text, type))
            text = 1e30
    return text

packages = ['PW', 'PP', 'NEB', 'BANDS', 'DOS', 'PROJWFC']
parameters = {}
defaults = {}
for package in packages:
    print("{0:=^60}".format(package))
    parameters[package] = {}
    defaults[package] = {}
    url = 'https://www.quantum-espresso.org/Doc/INPUT_%s.html'%package
    html_text = requests.get(url).text
    soup = BeautifulSoup(html_text, 'html.parser')
    sections = soup.find_all('table', attrs={ "border": "0"})
    for section in sections:
        section_header = section.find_all('h2')[0].text
        section_type = section_header.split(':')[0].strip()
        section_name = section_header.split(':')[1].strip().split()[0]
        print(section_type, section_name)
        tables = section.find_all('table', attrs={'width':"100%"})
        parameters[package][section_name] = {}
        defaults[package][section_name] = {}
        for table in tables:
            trs = table.find_all('tr')
            key = trs[0].find_all('th')[0].text.strip()
            if len(trs[0].find_all('th')) == 2:
                type = trs[0].find_all('th')[0].text.strip()
            else:
                type = trs[0].find_all('td')[0].text.strip()
            if len(trs) == 3:
                default = modify_text(trs[1].find_all('td')[1].text, type = type)
                options = trs[2].find_all('tt')
                option_list = []
                for option in options:
                    for text in option.text.split(','):
                        value = [modify_text(text, type = type)]
                        option_list = option_list + value
            else:
                default = modify_text('', type=type)
                options = trs[1].find_all('tt')
                option_list = []
                for option in options:
                    value = modify_text(option.text, type = type).split(',')
                    option_list = option_list + value
            parameters[package][section_name][key] = [type, option_list]
            defaults[package][section_name][key] = default

import pprint
pprint.pprint(parameters)
pprint.pprint(defaults)
with open('input_parameters.json', 'w') as f:
    json.dump(parameters, f, indent=2)
with open('input_parameters.json', 'w') as f:
    json.dump(default, f, indent=2)