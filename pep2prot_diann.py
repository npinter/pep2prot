import os
import sys
import time
import yaml
import pandas as pd
import numpy as np
import multiprocessing
from Bio import SeqIO
from termcolor import colored

global start_time


def print_welcome():
    os.system("color")
    print(colored(r"""
                 ____                  _   
 _ __   ___ _ __|___ \ _ __  _ __ ___ | |_ 
| '_ \ / _ \ '_ \ __) | '_ \| '__/ _ \| __|
| |_) |  __/ |_) / __/| |_) | | | (_) | |_ 
| .__/ \___| .__/_____| .__/|_|  \___/ \__|
|_|        |_|        |_|                               
DIA-NN edition""", "green"))
    print(colored("Niko Pinter - https://github.com/npinter/pep2prot", "green"))
    print(colored("v1.2", "green"))


def load_yaml(c_path):
    dir_files = os.listdir(c_path)
    yaml_list = []
    yaml_selected = None

    for file in dir_files:
        if file.endswith(".diann.yaml"):
            yaml_list.append(file)

    if len(yaml_list) == 1:
        yaml_selected = yaml_list[0]
        print("Using {} as input YAML file..".format(yaml_list[0]))
    elif len(yaml_list) > 1:
        yaml_loop = True
        while yaml_loop:
            print("Please select a YAML config file for pep2prot:")
            for yi, yaml_file in enumerate(yaml_list, start=1):
                print("  [{}] {}".format(yi, yaml_file))
            yaml_id_selected = input("Select ID: ")
            try:
                yaml_selected = yaml_list[int(yaml_id_selected) - 1]
                yaml_loop = False
            except ValueError:
                print(colored("Error: Type in a number!", "red"))
            except IndexError:
                print(colored("Error: ID not in list!", "red"))
    else:
        sys.exit(colored("Error: No YAML file found! Exit program!", "red"))

    return yaml_selected


def read_yaml(yaml_file):
    with open(yaml_file) as stream:
        return yaml.safe_load(stream)


def read_fasta(fasta_file):
    # create fasta data frame
    with open(fasta_file) as fasta:
        ids = []
        seqs = []
        name = []
        desc = []
        fasta_df_temp = pd.DataFrame()

        for seq_record in SeqIO.parse(fasta, "fasta"):
            ids.append(seq_record.id)
            seqs.append(seq_record.seq)
            name.append(seq_record.name)
            desc.append(seq_record.description)

        fasta_df_temp["Identifier"] = ids
        fasta_df_temp["Sequence"] = seqs
        fasta_df_temp["Name"] = name
        fasta_df_temp["Description"] = desc

    return fasta_df_temp


def annotate_peptides(jid, peptide_seqs, fasta_seqs, match_df_dict):
    match_peptide = []
    match_prot_id = []
    match_prot_name = []
    match_prot_description = []
    match_prot_seq = []
    match_prot_seq_start = []
    match_prot_seq_end = []
    match_prot_seq_aa_prv = []
    match_prot_seq_aa_last = []
    match_prot_seq_10_prv = []
    match_prot_seq_10_flw = []
    match_semi = []

    k = None
    len_peptide = len(peptide_seqs)
    proc_count = 0
    proc_count_perc = None

    # search peptide sequence in every fasta sequence
    for k, pep_seq in enumerate(peptide_seqs):
        if k > proc_count:
            if proc_count > 0:
                proc_count_perc = round(proc_count / (len_peptide / 100))
            else:
                proc_count_perc = 0

            print("Job #{}: {} of {} peptides ({}%)".format(jid, k, len_peptide, proc_count_perc))
            proc_count += int(len_peptide / 20)

        match_df = fasta_seqs[fasta_seqs["Sequence"].str.contains(pep_seq, regex=False)]

        if len(match_df):
            for m in range(0, len(match_df)):
                # write peptide sequence
                match_peptide.append(pep_seq)

                # write matched protein identifier
                ident_split = match_df["Identifier"].iloc[m].split("|")
                if len(ident_split) > 1:
                    match_prot_id.append(ident_split[1])
                else:
                    match_prot_id.append(ident_split[0])

                # write matched protein name
                prot_name = match_df["Identifier"].iloc[m].split("|")[-1]
                match_prot_name.append(prot_name)

                # write matched protein description
                try:
                    match_prot_description.append(
                        match_df["Description"].iloc[m].split(" ", maxsplit=1)[1].split("OS=")[0]
                    )
                except IndexError:
                    match_prot_description.append(
                        match_df["Description"].iloc[m].split(" ", maxsplit=1)[0]
                    )

                # write matched protein sequence
                prot_seq = match_df["Sequence"].iloc[m]
                match_prot_seq.append(prot_seq)

                # get peptide sequence start & end position in protein sequence
                prot_seq_start = match_df["Sequence"].iloc[m].index(pep_seq)
                prot_seq_end = prot_seq_start + len(pep_seq)
                match_prot_seq_start.append(prot_seq_start + 1)
                match_prot_seq_end.append(prot_seq_end)

                # get 1 and 10 previous and following amino acids of protein sequence
                prot_seq_10_prv = prot_seq[[prot_seq_start - 10
                                            if prot_seq_start > 10
                                            else None][0]:prot_seq_start]
                prot_seq_10_flw = prot_seq[prot_seq_end:[prot_seq_end + 10
                                                         if prot_seq_end + 10 <= len(prot_seq)
                                                         else None][0]]

                prot_seq_aa_prv = prot_seq_10_prv[-1:] if len(prot_seq_10_prv) > 1 else prot_seq_10_prv
                prot_seq_aa_last = pep_seq[-1:] if len(pep_seq) > 1 else pep_seq

                match_prot_seq_aa_prv.append(prot_seq_aa_prv)
                match_prot_seq_aa_last.append(prot_seq_aa_last)
                match_prot_seq_10_prv.append(prot_seq_10_prv)
                match_prot_seq_10_flw.append(prot_seq_10_flw)

                # semi-tryptic or not
                if prot_seq_aa_prv in ["K", "R"] or prot_seq_aa_last in ["K", "R"]:
                    if prot_seq_aa_prv in ["K", "R"] and prot_seq_aa_last not in ["K", "R"]:
                        match_semi.append("N")
                    elif prot_seq_aa_prv not in ["K", "R"] and prot_seq_aa_last in ["K", "R"]:
                        match_semi.append("C")
                    else:
                        match_semi.append("0")
                else:
                    match_semi.append("-1")

    # final job status output
    print("Job #{}: {} of {} peptides ({}%)".format(jid, k, len_peptide, proc_count_perc))

    # create match data frame
    match_df_temp = pd.DataFrame()
    match_df_temp["Peptide"] = match_peptide
    match_df_temp["Protein"] = match_prot_id
    match_df_temp["Name"] = match_prot_name
    match_df_temp["Description"] = match_prot_description
    match_df_temp["Sequence"] = match_prot_seq
    match_df_temp["Start"] = match_prot_seq_start
    match_df_temp["End"] = match_prot_seq_end
    match_df_temp["PreviousAA"] = match_prot_seq_aa_prv
    match_df_temp["LastAA"] = match_prot_seq_aa_last
    match_df_temp["Previous10"] = match_prot_seq_10_prv
    match_df_temp["Following10"] = match_prot_seq_10_flw
    match_df_temp["SemiOrNot"] = match_semi

    match_df_dict[jid] = match_df_temp


if __name__ == "__main__":
    start_time = time.time()
    manager = multiprocessing.Manager()
    manager_dict = manager.dict()

    # print welcome message
    print_welcome()

    # load input file paths from YAML
    config_dir_path = os.path.join(os.path.dirname(__file__), "config")
    pep2prot_yaml = read_yaml(os.path.join(config_dir_path, load_yaml(config_dir_path)))

    # peptide_tsv_path = pep2prot_yaml["peptide_tsv_path"]
    protein_expression_csv = pep2prot_yaml["protein_expression_csv"]

    fasta_db_path = pep2prot_yaml["fasta_db_path"]

    # set number of parallel jobs
    if pep2prot_yaml["num_jobs"] == 0:
        num_jobs = multiprocessing.cpu_count()
    else:
        num_jobs = pep2prot_yaml["num_jobs"]

    # create peptide and fasta dataframes
    peptide_df = pd.read_csv(protein_expression_csv, sep=";", header=0)
    fasta_df = read_fasta(fasta_db_path)

    # rename peptide column
    peptide_df = peptide_df.rename(columns={"{}".format(pep2prot_yaml["peptide_column_name"]): "Peptide"})

    # drop peptide duplicates
    peptide_df_unique = peptide_df.drop_duplicates(subset=["Peptide"])

    # setup all jobs
    match_df_multi = []
    peptide_df_max = 0
    peptide_df_step = int(np.ceil(len(peptide_df_unique) / num_jobs))

    # print job message
    print(colored("Mapping a total of {} peptides in {} jobs..".format(len(peptide_df_unique), num_jobs), "green"))

    jobs = []
    for i in range(0, num_jobs):
        peptide_df_min = peptide_df_max
        peptide_df_max += peptide_df_step

        if peptide_df_max > len(peptide_df_unique):
            peptide_df_max = len(peptide_df_unique)

        peptide_seqs_slice = peptide_df_unique["Peptide"][peptide_df_min:peptide_df_max]

        process = multiprocessing.Process(target=annotate_peptides,
                                          args=(i, peptide_seqs_slice, fasta_df, manager_dict))
        jobs.append(process)

    # start jobs
    for j in jobs:
        j.start()

    # validate if all jobs are finished
    for j in jobs:
        j.join()

    match_df_out = pd.concat(manager_dict.values(), ignore_index=True)

    # save dataframe to tsv file
    save_err = True
    while save_err:
        try:
            match_df_out.to_csv(os.path.join(pep2prot_yaml["output_folder"],
                                "pep2prot_mapped_peptides_{}.tsv".format(int(start_time))),
                                sep="\t",
                                index=False)
            save_err = False
        except PermissionError:
            print(colored("Error: Please close the \"pep2prot_mapped_peptides.tsv\" file!", "red"))
            input("Press \"ENTER\" to continue . . .")
        except FileNotFoundError:
            output_folder_name = os.path.basename(pep2prot_yaml["output_folder"])
            print(colored(
                "Error: Missing output folder! Creating {} in {}".format(output_folder_name,
                                                                         pep2prot_yaml["output_folder"]),
                "red"))
            os.mkdir(pep2prot_yaml["output_folder"])
            input("Press \"ENTER\" to continue . . .")

    end_time = time.time()
    print(colored("Runtime (total): {}".format(str(round((end_time - start_time) / 60, 2))), "green"))
