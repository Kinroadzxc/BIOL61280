import os

# Dictionary for ORF detection
search_pattern = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
                  "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S", "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
                  "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
                  "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
                  "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
                  "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
                  "TAA": ".", "TAC": "Y", "TAG": ".", "TAT": "Y", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
                  "TGA": ".", "TGC": "C", "TGG": "W", "TGT": "C", "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}
mutation_pattern = ('A', 'G', 'C', 'T')  # Possible mutations
CODON_LENGTH = 3  # the codon length is three


# Read in data and return the sequence as a string
def read_file(file_name):
    seq = ""  # initiate the sequence variable
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('>'):  # skip the line with name
                continue
            seq += line.strip()  # combine each line
    return seq


# Translate ORF
def translate(seq):
    if seq in search_pattern.keys():
        return search_pattern[seq]


# Mutate, translate and compare
def mutate_and_compare(seq, original_protein):
    mutation_list = [""]  # list containing all the single position mutation
    for i in range(len(seq)):
        if i < 3:  # ignore the start codon
            continue
        for nuc in mutation_pattern:
            tmp = ""
            if i % 3 == 2:
                tmp = seq[i-2:i] + nuc
            elif i % 3 == 1:
                tmp = seq[i-1] + nuc + seq[i+1]
            elif i % 3 == 0:
                tmp = nuc + seq[i+1:i+3]
            mutated = translate(tmp)  # cache the mutated sequence
            if mutated == '.':
                continue  # avoid stop codon
            # then compare the mutation with original protein sequence
            index = int(i / 3)
            if index >= len(original_protein):
                continue
            if original_protein[index] != mutated:
                result = [index + 1, original_protein[index], mutated]
                if result != mutation_list[-1]:  # avoid duplication
                    mutation_list.append(result)
    return mutation_list[1:]


# Write results into a fasta file
def output(out_list, gene_name):
    out_dir = "./output.fasta"  # the output path
    f_name = open(out_dir, 'a+')  # use 'a+' to avoid file overwrite
    for result in out_list:  # the format is "[gene name], [mutation index], [original aa], [mutated aa]"
        print(gene_name, result[0], result[1], result[2], file=f_name)


if __name__ == "__main__":
    p_list = []  # the list containing the protein names
    input_dir = "input3"

    # search for all the input folders, input file folders were named after the protein
    with os.scandir(input_dir) as entries:
        for entry in entries:
            # ignore system hidden files
            if str(entry.name).startswith('.'):
                continue
            # add all the folders to the protein list
            p_list.append(str(entry.name))

    # perform the mutation generator in each folder
    for path in p_list:
        # Gene sequence is waiting to be mutated while protein sequence is prepared for comparison
        gene = protein = ""  # initiate the sequence.
        with os.scandir(input_dir + '/' + path) as entries:  # read the sequences
            for entry in entries:
                if str(entry.name).startswith('.'):  # ignore system hidden files
                    continue
                if str(entry.name).startswith('seq'):  # nucleotide sequence is name like "seq_[geneName]"
                    gene = read_file(input_dir + '/' + path + "/" + str(entry.name))
                else:  # protein sequence is named like "P38398" or "O43823"
                    protein = read_file(input_dir + '/' + path + "/" + str(entry.name))
        if protein != "" and gene != "":  # perform the generator if two sequences have been loaded
            output(mutate_and_compare(gene, protein), path)
            print(path + " done")
