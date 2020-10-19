#Configuration dictionary for submitting jobs
import os

# default settings
xespressorc = {
    'queue':{
          # 'account': 'dcb',
          'time': '23:59:00',
        #   'mem-per-cpu': '5G',
          'nodes': '1',
          #'ntasks-per-node': '20',
          'cpus-per-task': '1',
          # 'constraint': 'mc',
          #'partition': 'empi',
          },
    'script': 'None',
      }

config_files = [os.path.join(os.environ['HOME'], '.xespressorc'),
            '.xespressorc']
for cf in config_files:
    if os.path.exists(cf):
        file = open(cf, 'r')
        xespressorc['script'] = file.read()
