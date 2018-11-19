#!/usr/bin/env python3
'''
CGE db checker
Checks the provided database folder (db's end with .fsa) and checks that there are no empty lines (removes if there
are), checks when split on delimiter that entries have the same number of items (warnings issued if not), and checks
for no duplicate entries in the database (warnings if not). Also does formatting according to standard in scikit.bio
for fasta output format.
input = database file to check (db's end with .fsa)
output = location of new file (new file will be checked for errors as well as modified for formatting)
delim = the expected delimiter that the database currently splits on, assumes all within folder use same delimiter
values = the number of expected values the delimiter will split into, assumes all within folder have same number
'''
import argparse
import os
import sys
import skbio

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Update pub mlst dbs')
    parser.add_argument("-in", "--input_database",
                        default="/srv/data/DB/cge_dbs/resfinder",
                        type=str)
    parser.add_argument("-out", "--out_database",
                        default="/srv/data/DB/cge_dbs/modified/resfinder",
                        type=str)
    parser.add_argument("-delim", "--delimiter",
                        default="_",
                        type=str)
    parser.add_argument("-values", "--expected_num_of_values",
                        default=2,
                        type=int)
    parser.add_argument("-file_type", "--file_extension",
                        default=".fsa",
                        type=str)

    args = parser.parse_args()
    return args

def check_cge_db(input_database, out_database, delimiter, expected_num_of_values, file_extension):
    sys.stdout.write("Check {} database location\n".format(input_database))
    if delimiter == "":
        delimiter = None
    if not os.path.isdir(out_database):
        try:
            os.makedirs(out_database)
        except Exception as e:
            sys.stderr.write("Error making output directory\n")
            return 1
    for database in sorted(os.listdir(input_database)):
        if database.endswith(file_extension):
            sys.stdout.write("Starting check on {}\n".format(database))
            with open(os.path.join(input_database, database), "r") as handle:
                fasta_map = {}
                number_of_errors = 0
                tmp = handle.read()
                if tmp.count("\n\n") > 0:
                    number_of_errors = tmp.count("\n\n")
                    sys.stdout.write("{} of empty lines are present, removing them\n".format(tmp.count("\n\n")))
                tmp = tmp.replace("\n\n", "\n").split("\n")  # Remove empty lines
                buffer = [i + "\n" for i in tmp]
                try:
                    entry_number = 0
                    fasta_db = {}
                    for seq in skbio.io.read(buffer, format="fasta"):
                        entry_number += 1
                        if len(seq.metadata['id'].split(delimiter)) != expected_num_of_values:
                            number_of_errors +=1
                            sys.stdout.write("Entry number {}, id:{}, expected {} items got {} items: [{}]\n".format(entry_number, seq.metadata['id'], expected_num_of_values, len(seq.metadata['id'].split(delimiter)), ','.join(seq.metadata['id'].split(delimiter))))
                        if seq.metadata['id'] in fasta_db:
                            number_of_errors +=1
                            sys.stderr.write("Duplicate header {} at entry {}\n".format(seq.metadata['id'], entry_number))
                        else:
                            fasta_db[seq.metadata['id']] = skbio.sequence.Sequence(seq.values.tostring().upper(), metadata={'id': seq.metadata['id']}) #  Note: This formats the length and ensures all letters are ALL CAPS                
                except Exception as e:
                    sys.stderr.write("Error with file, likely not fasta file. Error {}".format(str(e)))
                else:
                    try:
                        with open(os.path.join(out_database, database), "w") as output:
                            for entry in fasta_db: #  NOTE: python 3.6 dict order is guranteed on insertion order
                                fasta_db[entry].write(output)
                    except Exception as e:
                        sys.stderr.write("File creation error {}\n".format(str(e)))
                    else:
                        sys.stdout.write("Done creating output {} with {} errors. Formatting was ensured.\n\n".format(os.path.join(out_database, database), number_of_errors))
    return 0

if __name__ == "__main__":
    args = parse_args(sys.argv)
    check_cge_db(args.input_database, args.out_database, args.delimiter, args.expected_num_of_values, args.file_extension)
