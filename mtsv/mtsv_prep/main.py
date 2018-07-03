import subprocess
import sys
import argparse
import os
import datetime
import inspect
from ftplib import FTP
from time import sleep
import gzip
import tarfile
from multiprocessing import Pool, Queue, Process, Manager, RLock, freeze_support
from glob import iglob
from mtsv.commands import Command
from mtsv.parsing import make_sub_parser, get_global_config, ACTIONS, TYPES, add_default_arguments
from mtsv.mtsv_prep.MTSv_prune import *

from mtsv.utils import error, bin_path, specfile_path, specfile_read

from snakemake.workflow import Workflow, Rules, expand
import snakemake.workflow
from snakemake import shell

class Database(Command):
    config_section = ["DATABASE"]

    def __init__(self, params):
        print("running Database")
        super().__init__(params)
        print(self.params)


class CustomDB(Command):
    config_section = ["CUSTOM_DB"]

    def __init__(self, params):
        print("running custom database")
        super().__init__(params)
        print(self.params)

DEFAULT_DB_PATH = os.path.join(
    inspect.getfile(inspect.currentframe()),
    datetime.datetime.now().strftime("%b-%d-%Y"))

DEFAULT_PARTITIONS = ["2,2157", "10239,12884", "28384",
                      "2759-33090,4751,7742", "33090",
                      "4751", "7742-9443,9397,9913,9615",
                      "9443-9606", "9397", "9913", "9615", "9606"]

COMMANDS = {
    "database" : Database,
    "custom_db": CustomDB,
    }

FILE_LOCK = RLock()

def oneclickdl(args):
    fin = set()
    for db in ["genbank", "Complete_Genome", "Chromosome", "Scaffold"]:
        if args.path and os.path.isfile(os.path.join(args.path, "artifacts","{0}.fas")):
           fin.add(db)
    dbs = set(args.includedb).difference(fin)
    args.includedb = list(dbs)

    return pull(
        path=args.path,
        thread_count=args.threads,
        databases=args.includedb)

def decompression(x, path):
    output = subprocess.PIPE
    if not os.path.isfile(x.strip(".gz")):
        subprocess.run(['gunzip',x], stderr=output, stdout=output)
        if output:
            with FILE_LOCK:
                with open(os.path.join(path, "artifacts/decompression.log"), "a") as out:
                    out.write("{0}\n".format(x))
    return x.strip(".gz")

def oneclickbuild(args):
    with open(os.path.join(args.path, "artifacts/decompression.log"), "w" ):
        pass
    pool = Pool(args.threads)
    pool.starmap(decompression, [(os.path.abspath(x),
                                  args.path) for x in iglob(os.path.join(args.path,"flat_files/*.gz"))])
    for fp in iglob(os.path.join(args.path,"artifacts/","*_ff.txt")):
        db = list(os.path.split(fp))
        db[1] = db[1].strip().replace("_ff.txt",".fas")
        db = os.path.abspath("{0}{1}{2}".format(db[0],os.sep,db[1]))
        with open(os.path.abspath(fp), "r" ) as file:
            temp = file.readlines()

        with open(os.path.abspath(fp), "w") as file:
            for x in temp:
                if x.strip().rsplit(".",1)[1] == "gz":
                    x= x.strip().rsplit(".",1)[0]
                file.write("{0}\n".format(x))

        if os.path.isfile(db):
            continue
        build_db(os.path.abspath(fp), db, os.devnull, os.devnull, args.threads, os.devnull)

    arguments = oneclickjson(args.path)

    pool.starmap(acc_serialization, [(argument['acc-to-taxid-paths'], argument['fasta-path'],
                                      argument['taxdump-path']) for argument in arguments ])

    # shutil.rmtree(os.path.join(args.path, "flat_files" ))

def mapper(x):
    return clip(*x)

def partition(args):
    partition_list = set()
    for db in args.customdb:
        db = db.strip().lower().replace(" ","_")

        arguments = parse_json(os.path.join(args.path, "artifacts/{0}.json".format(db)))
        for prt in args.partitions:
            try:
                temp = prt.split("-")
                if len(temp) == 2:
                    inc = set(temp[0].split(","))
                    exc = set(temp[1].split(","))
                    path = os.path.join(args.path, "indices", db, "{0}-{1}".format("_".join(sorted(inc)), "_".join(sorted(exc)) ))
                    os.makedirs(path, exist_ok=True)
                    prt = "{0}-{1}".format("_".join(sorted(inc)), "_".join(sorted(exc)))
                else:
                    inc = set(temp[0].split(","))
                    exc = set()
                    path = os.path.join(args.path, "indices", db, "{0}".format("_".join(sorted(inc)) ))
                    os.makedirs(path, exist_ok=True)
                    prt = "{0}".format("_".join(sorted(inc)) )
                partition_list.add( ( list(inc), args.rollup_rank, list(exc), os.path.join(path,
                                        "{0}.fas".format(prt)),arguments['minimum-length'],
                                         arguments['maximum-length'], arguments["fasta-path"],
                                         arguments["serialization-path"], args.debug  ) )
            except OSError:
                print("Partion folder {0} exists please use --overwrite to repartition".format(path))

    p = Pool(args.threads)
    ret_list = p.starmap(clip, partition_list)

    return ret_list, args

def chunk(file_list, args):
    dir_set = set()
    for fp in file_list:
        if not os.path.isfile("_0.".join(fp.rsplit(".", 1))):
            subprocess.run("{2} --input {0} --output {1} --gb 2".format(fp, os.path.dirname(fp), bin_path('mtsv-chunk')).split() )
        dir_set.add(os.path.dirname(fp))
    return list(dir_set)

def snake(args):
    chunk_path = bin_path('mtsv-chunk')
    fm_build_path = bin_path('mtsv-build')
    print(chunk_path)
    print(fm_build_path)
    print(args)

    print(partition(args))
    # for i in args:
    #     print(i)
    workflow = Workflow(
        "__file__",
        overwrite_workdir=args.working_dir)
    if args.cluster_cfg is not None:
        workflow.cluster_cfg = args.cluster_cfg
    snakemake.workflow.rules = Rules()
    snakemake.workflow.config = dict()




    # @workflow.rule(name='chunking')
    # @workflow.docstring("""Fasta chunking""")
    # @workflow.input(os.path.join("test", "{}.index"))
    # @workflow.output(os.path.join(cmd.params['binning_outpath'], "{index}.bn"))
    # @workflow.params(call=chunk_path, args=args.chunk_size)
    # @workflow.message("Executing Fasta chunking on the follow files {input}.")
    # @workflow.log(os.path.join(args.path,"artifacts/logs/""{index}.log"))
    # @workflow.threads(1)
    # @workflow.shellcmd("{params.call} --input {input} --output{output} --gb {params.args} > {log} 2>&1")
    # @workflow.run
    # def __rule_chunking(input, output, params, wildcards,
    #                    threads, resources, log, version, rule,
    #                    conda_env, singularity_img, singularity_args,
    #                    use_singularity, bench_record, jobid, is_shell):
    #     shell(
    #         "{params.call} --input {input} --output{output} --gb {params.args} > {log} 2>&1"
    #     )
    #
    # workflow.check()
    # print("Dry run first ...")
    # workflow.execute(dryrun=True, updated_files=[])
    # return workflow

def fm_build(dir_list):
    fm_list = []
    for directory in dir_list:
        for fp in iglob(os.path.join(directory, "*.fasta")):
            out_file = os.path.join(directory, "{0}.index".format(os.path.basename(fp).split(".")[0]))
            subprocess.run("{2} --fasta {0} --index {1}".format(os.path.abspath(fp), out_file, bin_path('mtsv-build')).split() )
            fm_list.append(out_file)
    return fm_list

def oneclickfmbuild(args, is_default):
    to_link = fm_build(chunk(partition(args)))


def json_updater(args):

    for json_path in iglob(os.path.join(args.path, "artifacts/*.json")):
        base = os.path.basename(json_path)
        if base == "exclude.json":
            continue
        base = base.rsplit(".", 1)[0]

        with open(json_path, "r") as file:
            params = json.load(file)
        fm_all = []
        params['fm-paths'] = {}
        for path in iglob(os.path.join(args.path, "indices", base, "**/*.index")):
            folder = os.path.basename(os.path.dirname(path))
            try:
                params['fm-paths'][folder].append(os.path.abspath(path))
            except KeyError:
                params['fm-paths'][folder] = [os.path.abspath(path)]
            fm_all.append(os.path.abspath(path))
        params['fm-index-path'] = fm_all
        params['partition-path'] = []

        for path in iglob(os.path.join(args.path, "fastas", base, "*.fas")):
            # folder = os.path.basename(os.path.dirname(path))
            # print(path)
            # try:
            #     params['partition-path'][folder]
            # except KeyError:
            #     params['partition-path'][folder] = []
            params['partition-path'].append(os.path.abspath(path))

        params['mtsv-tree'] = os.path.abspath(os.path.join(args.path, "artifacts","tree.index"))
        with open(json_path, "w") as file:
            json.dump(params, file, sort_keys=True, indent=4)

def make_json_rel(args):
    rm_path = os.path.abspath(args.path)
    for name in ["genbank", "Complete_Genome","Chromosome","Scaffold"]:
        try:
            arguments = parse_json(os.path.join(args.path, "artifacts","{0}.json".format(name)))
            try:
                for i, abs_path in enumerate(arguments['acc-to-taxid-paths']):
                    arguments['acc-to-taxid-paths'][i] = os.path.relpath(abs_path, rm_path)
            except KeyError:
                pass
            try:
                for j in arguments['fm-paths'].keys():
                    for i, abs_path in enumerate(arguments['fm-paths'][j]):
                        arguments['fm-paths'][j][i] = os.path.relpath(abs_path, rm_path)
            except KeyError:
                pass
            try:
                for j in arguments['partition-path'].keys():
                    for abs_path in arguments['partition-path'][j]:
                        arguments['partition-path'][j] = os.path.relpath(abs_path, rm_path)
            except KeyError:
                pass

            for key in arguments.keys():
                try:
                    if rm_path in arguments[key]:
                        arguments[key] = os.path.relpath(arguments[key], rm_path )
                    # else:
                except TypeError:
                    continue
        except FileNotFoundError:
            continue
        with open(os.path.join(args.path, "artifacts","{0}.json".format(name)), "w") as file:
            json.dump(arguments, file, sort_keys=True, indent=4)



def make_json_abs(args):
    for name in ["genbank", "Complete_Genome","Chromosome","Scaffold"]:
        try:
            arguments = parse_json(os.path.join(args.path, "artifacts","{0}.json".format(name)))
            try:
                for i, abs_path in enumerate(arguments['acc-to-taxid-paths']):
                    arguments['acc-to-taxid-paths'][i] = os.path.abspath(os.path.join(args.path,abs_path))
            except KeyError:
                pass
            try:
                # keys = list(arguments['fm-paths'].keys())
                for j in arguments['fm-paths'].keys():
                    for i, abs_path  in enumerate(arguments['fm-paths'][j]):
                        arguments['fm-paths'][j][i] = os.path.abspath(os.path.join(args.path,abs_path))
            except KeyError:
                pass
            try:
                # keys = list()
                for j, val in enumerate(arguments['partition-path']):
                    # for abs_path in arguments['partition-path'][j]:
                    arguments['partition-path'][j] = os.path.abspath(val)
            except KeyError:
                pass

            for key in arguments.keys():
                try:
                    if os.path.isfile(os.path.abspath(os.path.join(args.path,arguments[key]))):
                        arguments[key] = os.path.abspath(os.path.join(args.path, arguments[key]))
                except TypeError:
                    continue
        except FileNotFoundError:
            continue
        with open(os.path.join(args.path, "artifacts","{0}.json".format(name)), "w") as file:
            json.dump(arguments, file, sort_keys=True, indent=4)


def setup_and_run(parser):
    if sys.argv[1] == "json_update":
        args = parser.parse_known_args()[0]
        json_updater(args)
        make_json_abs(args)

    else:
        args = parser.parse_known_args()[0]

        try:

            if args.cmd_class == Database:
                if args.download_only:
                    args.path = os.path.abspath(oneclickdl(args))

                elif args.build_only:
                    if args.path and  os.path.isdir(args.path):
                        oneclickbuild(args)
                    else:
                        print("A valid path was not specified")
                else:
                    args.path = os.path.abspath(oneclickdl(args))
                    oneclickbuild(args)

            elif args.cmd_class == CustomDB:
                oneclickfmbuild(args, args.partitions == DEFAULT_PARTITIONS)
                make_json_abs(args)

            # try:
            #     make_json_rel(args)
            # except:
            #     pass

        except AttributeError:
            sys.argv[1] = "database"
            args = parser.parse_known_args()[0]
            args.path = os.path.abspath(oneclickdl(args))
            oneclickbuild(args)
            path = args.path
            sys.argv[1] = "custom_db"
            args = parser.parse_known_args()[0]
            args.path = path
            oneclickfmbuild(args, args.partitions == DEFAULT_PARTITIONS)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    

    parser = argparse.ArgumentParser(
        prog="mtsv-setup",
        description="Download and build sequence databases and indices",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(
        title="commands", metavar="COMMAND",
        help="Setup Commands"
    )

    for command, cmd_class in COMMANDS.items():
        make_sub_parser(
            subparsers, command, cmd_class
        )
    p = subparsers.add_parser("json_update")
    p.add_argument("--path")
    p = subparsers.add_parser("oneclick")
    for command, cmd_class in COMMANDS.items():
        for arg, desc in get_global_config(cmd_class.config_section).items():
            if "_meta" in arg:
                continue
            if 'type' in desc:
                desc['type'] = TYPES[desc['type']]
            if 'default' in desc and 'help' in desc:
                desc['help'] += " (default: {})".format(desc['default'])
            if 'action' in desc and desc['action'] in ACTIONS:
                desc['action'] = getattr(
                    sys.modules[__name__], desc['action'])
            arg = "--{}".format(arg)
            if 'positional' in desc:
                del desc['positional']
            try:
                p.add_argument(
                    arg, **desc
                )
            except argparse.ArgumentError:
                continue
        try:
            add_default_arguments(p)
        except argparse.ArgumentError:
            pass

    if len(argv)==1:
        parser.print_help(sys.stdout)
        sys.exit(1)
    try:
        setup_and_run(parser)
    except KeyboardInterrupt:
        error("\n-- Stopped by user --", exception=False)


if __name__ == "__main__":
    main()
