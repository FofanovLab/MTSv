import subprocess
from io import BytesIO
import argparse
import fnmatch
import os, datetime
from ftplib import FTP
from time import sleep
import gzip
import tarfile
from multiprocessing import Pool, Queue, Process, Manager, RLock
import pickle, json


def cmd(x):
    with open(os.devnull, "w") as null:
        subprocess.run(x.split(), stderr=null, stdout=null)
    return x
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TaxClipper is intended to be used to parse sequences based on NCBI taxid")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-oc", "--oneclick", "-oneclick",action='store_true')
    group.add_argument("-c", "--custom", "-custom", action='store_true')


    group = parser.add_argument_group()
    group.add_argument("-p","--path","-path", nargs=1,
                       help="Path to dated folder containing artifacts")

    group.add_argument("-prt","--partitions","-partition", nargs='*',
                       help="Space seperated list of '-' delimited tuples of comma seperated TaxIDs for inclusion and exclusion\n"
                            "ex. 1386-1392,1396 would create a FM-index of the genus Bacillus without B. anthracis and B. cereus\n"
                            "if no exlusion desire no '-' is needed")
    group.add_argument("-min","--minimum-length","-minimum-length", type=int,
                       help="Integer for minimum length of sequences to include")
    group.add_argument("-max","--maximum-length","-maximum-length", type=int,
                       help="Integer for maximum length of sequences to include")
    group.add_argument("-rur","--rollup-rank","-rollup-rank",
                       help="NCBI rank to set sequence headers i.e. None, genus, family et cetera\n (default: species)")

    group.add_argument("-t", "--threads", "-threads", type=int,
                       help="Specify total threads to spawn in FM-index creation")
    group.add_argument("-ow", "--overwrite", "-overwrite", help="Specify total threads to spawn in FM-index creation",
                       action='store_true')

    group.add_argument("-cdb", "--customdb", "-customdb", nargs='*',
                       help="Space delimited list of the following NCBI source: Chromosome, Scaffold, genbank. By default only sequences from")

    args = parser.parse_args()
    if args.threads:
        threads = args.threads
    else:
        threads = 1

    if args.oneclick:

        if args.path:
            path = args.path
        else:
            path = datetime.datetime.now().strftime("%b-%d-%Y")
            with open(os.devnull, "w") as null:
                subprocess.run("python MTSv_prune.py -oc -t {0}".format(threads).split(), stdout=null,stderr=null)
        partitions = []
        if args.customdb:
            dbs = args.customdb
        else:
            dbs = ["genbank"]
        if args.partitions:
            default_parts = args.partitions
        else:
            #Chunk 1:All Bacteria, Chunk 2 All Viruses Chunke 3 Other Sequence Chunk 4 Eukaryotes minus plants,fungi vertebrates
            #Chunk 5 Plants Chunk 6 Fungi Chunk 7 Vertebrates minus primates bats cow dog Chunk 8 Primates minus human
            #Chunk 9 Bats Chunk 10 Cows Chunk 11 Dogs Chunk 12 Humans
            default_parts = ["2,2157", "10239,12884", "28384", "2759-33090,4751,7742",
                             "33090", "4751", "7742-9443,9397,9913,9615", "9443-9606"
                             "9397", "9913", "9615","9606"]
        fasta_list = []
        for db in dbs:
            for part in default_parts:
                out_folder = part.replace(",","_")
                try:
                    os.makedirs(os.path.join(path,"indices",db,out_folder), exist_ok=args.overwrite)
                    temp = part.split("-")
                    if len(temp) == 2:
                        partitions.append("python MTSv_prune.py -cp {0}.p  -o {1}.fasta -txi {2} -txe {3}".format(
                            os.path.join(path,"artifacts",db), os.path.join(path,"indices",db,out_folder, db),
                            " ".join(temp[0].split(",")), " ".join(temp[1].split(","))))
                    elif len(temp) == 1:
                        partitions.append("python MTSv_prune.py -cp {0}.p -o {1}.fasta -txi {2}".format(
                            os.path.join(path, "artifacts", db), os.path.join(path, "indices", db, out_folder, db),
                            " ".join(temp[0].split(","))))
                    else:
                        raise ValueError
                    fasta_list.append("{0}.fasta".format(os.path.join(path,"indices",db,out_folder, db)))
                except OSError:
                    print("FM-Index Directory Appears to Exist Remove {0} or use --overwrite flag if a rebuild is required.".format(os.path.abspath(os.path.join(path,"indices",db,out_folder))))
                    continue
        with Pool(threads) as p:
            ret = p.map(cmd, partitions)
        for i in ret:
            print(i)
