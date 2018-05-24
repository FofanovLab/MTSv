import json
import yaml
from collections import OrderedDict

def from_json(json_str):
	"""Reads a JSON string into an OrderedDict.

	:param json_str: a JSON string
	:returns: an OrderedDict of the JSON contents
	"""
	argsdict = json.loads(json_str, object_pairs_hook=OrderedDict)
	return argsdict

def from_yaml(yaml_str):
	"""Reads in a string of YAML into an OrderedDict.

	:param yaml_str: the contents of a YAML file
	:returns: an OrderedDict of the YAML contents
	"""
	def _dict_constructor(loader, node):
		return OrderedDict(loader.construct_pairs(node))
	_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
	yaml.add_constructor(_mapping_tag, _dict_constructor)

	argsdict = yaml.load(yaml_str)
	return argsdict

