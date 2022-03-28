import argparse
from ruamel import yaml
from ruamel.yaml.scalarstring import PreservedScalarString as pss
import sys
import pdb


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-cwl', action='store', dest='cwl', help='Input cwl file',required=True)
parser.add_argument('-r', '--readme', action='store', dest='readme', help='Readme file to insert into workflow/tool doc, if applicable', required=False)
parser.add_argument('-n', '--id-name', action='store', dest='id_name', help='Short app ID link name to use, i.e. kfdrc-align-wf', required=False)
parser.add_argument('-l', '--label', action='store', dest='label', help='User-friendly label to add to tool/workflow cwl, if needed', required=False)
parser.add_argument('-t', '--tags', action='store', dest='tags', help='Seven bridges tags file, as csv string, ex RNASEQ,FUSION', required=False)
parser.add_argument('-f', '--files', action='store', dest='files', help='Cavatica-style tsv manifest with file ID, file name, associated cwl input key, and, optionally, secondaryFile names and IDs', required=False)
parser.add_argument('-p', '--publisher', action='store', dest='pub', help='Publisher name', required=False, default="KFDRC")
args = parser.parse_args()


def update_file_paths(manifest, yaml_obj):
    mf = open(manifest)
    head = next(mf) # fields: fid,name,key,sfids,snames
    # see if more than one file exists for input keys - will determine if values is to be File or File[]
    in_dict = {}
    for entry in mf:
        info = [None] * 5
        data = entry.rstrip('\n').split('\t')
        info[:len(data)] = data
        fid,name,key,sfids,snames = info
        if sfids in [None,"None"]:
            dat = {"class":"File","path":fid,"name":name}
        else:
            assert len(sfids.split(',')) == len(snames.split(',')), "Lengths of secondaryFile IDs and names MUST BE equal"
            sfidsList = sfids.split(',')
            snamesList = snames.split(',')
            sfileList = []
            for i in range(len(sfidsList)):
                sfileList.append({"class":"File","path":sfidsList[i],"name":snamesList[i]})
            dat = {"class":"File","path":fid,"name":name,"secondaryFiles":sfileList}
        if key not in in_dict:
            in_dict[key] = []
        in_dict[key].append(dat)
        for key in in_dict:
            if len(in_dict[key]) == 1:
                yaml_obj['inputs'][key]["sbg:suggestedValue"] = in_dict[key][0]
            else:
                yaml_obj['inputs'][key]["sbg:suggestedValue"] = in_dict[key]

# round tripper preservers order and formatting of keys and values
data = yaml.load(open(args.cwl), yaml.RoundTripLoader, preserve_quotes=True)
if args.files:
    update_file_paths(args.files, data)
# check for license, pub
if 'sbg:license' not in data:
    data['sbg:license'] = "Apache License 2.0"
if 'sbg:publisher' not in data:
    data['sbg:publisher'] = args.pub
if args.tags:
    data['sbg:categories'] = []
    for tag in sorted(args.tags.split(',')):
        data['sbg:categories'].append(tag)
key_list = list(data.keys())
if args.label:
    data.insert(key_list.index('id')+1, 'label', args.label)
    data['label'] = args.label

if args.readme:
    rm = open(args.readme)
    rm_str = rm.read()
    rm.close()
    # code from should improve readme appearance
    pss_rm_str = pss(rm_str)
    key_list = list(data.keys())
    if 'doc' not in data:
        try:
            data.insert(key_list.index('label')+1, 'doc', pss_rm_str)
        except Exception as e:
            sys.stderr.write(str(e) + "\nFailed to add doc field after label field, trying after id\n")
            data.insert(key_list.index('id')+1, 'doc', pss_rm_str)
    else:
        data['doc'] = pss_rm_str
if args.id_name:
    key_list = list(data.keys())
    if 'id' not in data:
        data.insert(key_list.index('class')+1, 'id', args.id_name)
    else:
        data['id'] = args.id_name

yaml.dump(data, sys.stdout, Dumper=yaml.RoundTripDumper, default_flow_style=False)
