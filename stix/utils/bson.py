import json
import numpy 
class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(CustomEncoder, self).default(obj)
        
def dict_to_json(data_dict):
    encoded = json.dumps(data_dict,cls=CustomEncoder)
    return  json.loads(encoded)
