import os #os is for listing files/maps in a directory


#In my case all files are csv files (separated with ;), except for the vcf's.

which_genes_are_effectors_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/snpeff_fasta/effectors_and_crns.csv"
#For later annotation if you use all genes or all genes with a signal peptide. If you've already predicted effectors
#Could also be used for any other list of genes with a certain annotation.
genes_you_are_interested_in_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/snpeff_fasta/all_genes.csv"
#These are all genes for which you will look if there is a pattern.
map_with_your_vcfs_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/snpeff_fasta/homo_vcfs"
#These are the vcf files, already annotated (in genes) and selected for homozygosity.
output_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/snpeff_fasta/summary_tables/"
#Directory in which you output will be placed
resistance_pattern_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/snpeff_fasta/summary_tables/rpf_table.csv"
#A file with information about strains of your crop and which isolates they are resistant to
presence_absence_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/internship_sjoerd/pres_abs_correctNames.csv"
#If you have this, a file with presence and absence of genes
extra_annotation_dir = "C:/Users/sjoerd/Documents/MCLS/Major internship/snpeff_fasta/Pfam_annotations.csv"
#In my case Pfam predicted domains.

list_of_isolates = ["Pfs1","Pfs2","Pfs3","Pfs4","Pfs5","Pfs6","Pfs7","Pfs8","Pfs9","Pfs10","Pfs11","Pfs12","Pfs13","Pfs14","Pfs15","Pfs16","A0340","ES1314","F-0543","NL0539","US1135","US1322","US1324","US1509"]
#List of all isolates of which you have a vcf, in order of the output (this was easier to just type out for me)
pfs_only = ["Pfs1","Pfs2","Pfs3","Pfs4","Pfs5","Pfs6","Pfs7","Pfs8","Pfs9","Pfs10","Pfs11","Pfs12","Pfs13","Pfs14","Pfs15","Pfs16"]
#List of only the isolates of which you know which resistances they can break, for pattern recognition. In order


#A. 
#This field is optional. It makes a list of all effectors (in hashname). Output: effector_list = ["hashname1";"hashname2"; etc.] 
effectors_file = open(which_genes_are_effectors_dir,"r")
effector_list = []
for gene_information in effectors_file: 
    gene_information = gene_information.split(";") 
    if gene_information[2] != '': 
        effector_list.append(gene_information[2])
effectors_file.close()

#B.
peptide_changes_dict = {} 
#peptide_changes_dict = {gene of interest : {peptide change:[isolates]}}
file_with_genes_of_interest = open(genes_you_are_interested_in_dir,"r")    
gene_name_hash_dict = {} #Output: {"pfs1|...":"hashname1"}
genes_of_interest = [] #Just a list with all genes
goi_dict = {} # All relevant information from the genes of interest file in a dict, key = hashname, value = gene_name;length;hashname;SignalP;RxLR;EER;WY;etc;
for line in file_with_genes_of_interest:
    line = line.split(";")
    if line[2]!= '':
        gene_name_hash_dict[line[0]]=line[2] #So for all genes, we now link the pfs1|... name to the hashname
        genes_of_interest.append(line[2])
        goi_dict[line[2]] = str(line[0]+";"+line[1]+";"+line[2]+";"+line[7]+";"+line[16]+";"+line[18]+";"+line[12]+";"+line[13]+";"+line[14]+";"+line[21]+";")
        peptide_changes_dict[line[2]]={}
		
#!!!!!!!!!!!!!! This bit below is dependend on A. Can be silenced. Needed for the summary table (is this an effector)
        if line[2] != " Hashname":
            if line[2] in effector_list:
                goi_dict[line[2]]+= "yes;"
            else:
                goi_dict[line[2]]+= "no;"
effectors_file.close()


list_of_homo_files = os.listdir(map_with_your_vcfs_dir) 
# This makes a list of the names of all files in that directory

#C.
# Dependend on B.
#This loops to every vcf file and finds which amino acids in which peptides are changed
#It looks for every line: Is there a gene of interest in there, is there an amino acid change
#Output location_change_dict ={hashname_Gene:{peptide_change:[Isolate;possibly_more_isolates]}}
for isolate in list_of_homo_files:
    location_change_dict = {} 
    isolate_file = open(map_with_your_vcfs_dir+isolate,"r")
    for line in isolate_file:
        if line[0] != "#":
            for gene in genes_of_interest:
                if gene in line:
                    isolate_information = line.split("|")
                    for piece_of_information in isolate_information:
                        if "p." in piece_of_information: #so if there is a peptide change
                            if piece_of_information[2:] in peptide_changes_dict[gene]: #so everything except for p.
                                peptide_changes_dict[gene][piece_of_information[2:]].append(isolate.strip("_MOHIHOMO.ann.vcf"))
                            else:
                                peptide_changes_dict[gene][piece_of_information[2:]] = []
                                peptide_changes_dict[gene][piece_of_information[2:]].append(isolate.strip("_MOHIHOMO.ann.vcf"))
    isolate_file.close()
    
#D.
#Dependend on B. & C.
#So this is where we make an output file with for every single SNP the pattern in which it is found.
annotated_table = open(output_dir+"all_signalp_annotated.csv", "w")
header = ";;" #This will become the header of the annotated table, starting with two open cells (excel). Could also go for "Gene name;Peptide change;"
mutations_dict = {} #Mutations_dict={Gene_mutation:[-;-;x;x;x;-;x;-;x;]} (=pattern of this specific snp/mutation)
for isolate in list_of_isolates:
    header+= isolate + ";"
for gene in peptide_changes_dict:  #{gene of interest : {peptide change:[isolates]}}
#    if peptide_changes_dict[gene]!={}:
#        annotated_table.write(header+"\n") #If the gene has any SNPs, start new line, else: Skip
    for mutation in peptide_changes_dict[gene]:
        write_line = gene+";"+mutation+";" #First two columns of the table are gene name and aa change
        mutations_dict[gene+"_"+mutation]= '' 
        for isolate in list_of_isolates: #for every isolate in order we look which snps are there, resulting in a pattern (x==snp, -==same as reference)
            if isolate in peptide_changes_dict[gene][mutation]:
                write_line+="x;"
                if isolate in pfs_only: #So if the isolate is not a US etc. version
                    if str(gene+"_"+mutation) not in mutations_dict: #if we haven't found this snp earlier, start this dict
                        mutations_dict[str(gene+"_"+mutation)]="x;"
                    else:
                        mutations_dict[str(gene+"_"+mutation)]+="x;" #if we have already found this variant, put x
            else:
                write_line += "-;" #There is a variant other than the reference, but it is not found in this isolate.
                if isolate in pfs_only:
                    if str(gene+"_"+mutation) not in mutations_dict:
                        mutations_dict[str(gene+"_"+mutation)]="-;" #We should ofcourse never find this
                        print("What's happening here??")
                    else:
                        mutations_dict[str(gene+"_"+mutation)]+="-;"
        if "x" in write_line:#We only put in genes that have at least one variant that is different from the reference
            annotated_table.write(write_line+"\n")
annotated_table.close()

#E.
# Dependend on the availability of a resistance table.
#Here we make a similar pattern with the resisntances as we have done with the variants, so they can easily be compared.
rpf_table = open(resistance_pattern_dir,"r")
patterns = []
# Just a list with all patterns that are in this table
pattern_dict = {}
#A dict with the pattern as a key and the spinach line as value
for row in rpf_table:
    pattern = row[2:]
    pattern = pattern.strip()
    pattern+= ";"
    pattern_dict[pattern] = row[0]
    patterns.append(pattern)
rpf_table.close()


#F.
# Dependend on B. C. D.
#Simplified, Finding if there are homozygote mutations in the same gene that follow a pattern
simplified_mutation_dict = {}
for gene in peptide_changes_dict:
    gene_isolate_list = []
    for mutation in peptide_changes_dict[gene]:
        for isolate in peptide_changes_dict[gene][mutation]:
            gene_isolate_list.append(isolate)
    for pfs in pfs_only:
        if pfs in gene_isolate_list:
            if gene not in simplified_mutation_dict:
                simplified_mutation_dict[gene] = "x;"
            else:
                simplified_mutation_dict[gene] += "x;"
        else:
            line += "-;"
            if gene not in simplified_mutation_dict:
                simplified_mutation_dict[gene] = "-;"
            else:
                simplified_mutation_dict[gene] += "-;"

#G.
# Dependend on E.
#Here I compare the presence absence with the patterns
absence_pres_file = open(presence_absence_dir,"r")
pres_abs_dict = {}
for gene in absence_pres_file:
    line = []
    gene = gene.split(";")
    for value in gene:
        if value == "1":
            line.append("-")
        elif value == "0":
            line.append("x")
        else:
            line.append(value)
    if "x" in line:
        pres_abs_dict[line[0]]=line[5]+";"+line[13]+";"+line[14]+";"+line[15]+";"+line[16]+";"+line[17]+";"+line[18]+";"+line[19]+";"+line[20]+";"+line[6]+";"+line[7]+";"+line[8]+";"+line[9]+";"+line[10]+";"+line[11]+";"+line[12]+";"
absence_pres_file.close()

#H.
# Depended on B. C. D. E. F. G.
#Simplified mutations together with presence_absence
simple_pres_mutation_dict = {}
for gene in simplified_mutation_dict:
    if "x" in simplified_mutation_dict[gene]:
        if gene in pres_abs_dict:
            if "x" in pres_abs_dict[gene]:
                simple_pres_mutation_dict[gene] = ""
                for r in range(0,len(simplified_mutation_dict[gene])):
                    r = int(r)
                    if simplified_mutation_dict[gene][r] == pres_abs_dict[gene][r]:
                        simple_pres_mutation_dict[gene]+=simplified_mutation_dict[gene][r]
                    else:
                        if simplified_mutation_dict[gene][r] == "x":
                            simple_pres_mutation_dict[gene]+=("x")
                        elif pres_abs_dict[gene][r] == "x":
                            simple_pres_mutation_dict[gene]+=("x")
                        else:
                            simple_pres_mutation_dict[gene]+=("-")


#I.
# Depended on B. C. D. E. F. G.
#Mutations together with presence_absence
pres_mutation_dict = {}
for mutation in mutations_dict:
    basename = gene.strip("_")
    if "x" in mutations_dict[mutation]:
        if basename in pres_abs_dict:
            if "x" in pres_abs_dict[basename]:
                pres_mutation_dict[mutation] = ""
                for r in range(0,len(mutations_dict[mutation])):
                    r = int(r)
                    if mutations_dict[mutation][r] == pres_abs_dict[basename][r]:
                        pres_mutation_dict[mutation]+=(mutations_dict[mutation][r])
                    else:
                        if mutations_dict[mutation][r] == "x":
                            pres_mutation_dict[mutation]+=("x")
                        elif pres_abs_dict[basename][r] == "x":
                            pres_mutation_dict[mutation]+=("x")
                        else:
                            pres_mutation_dict[mutation]+=("-")

#J.
# Depending on B
#Adding Pfam domains
pfam_file = open(extra_annotation_dir, "r")
pfam_dict = {}


for annotation_line in pfam_file:
    annotation_line = annotation_line.split(";")
    gene_name = annotation_line[0].replace("_","|")
    if gene_name in gene_name_hash_dict:
        if gene_name_hash_dict[gene_name] not in pfam_dict:
            pfam_dict[gene_name_hash_dict[gene_name]] = str(annotation_line[4]+";"+annotation_line[5])
        else:
            pfam_dict[gene_name_hash_dict[gene_name]] += str(";"+annotation_line[4]+";"+annotation_line[5])

for gene in genes_of_interest:
    if gene not in pfam_dict:
        pfam_dict[gene] = ""


#K.
#This will make a final summary of everthing above.If you don't use statemant A. or J., the header should be adjusted.
final_summary_file = open(output_dir + "final_summary.csv","w")
final_summary_file.write(str(goi_dict[" Hashname"]+"Predicted effector;Amino acid change;RPF;What kind of pattern;pfam domain;pfam annotation;pfam domain;pfam annotation;pfam domain;pfam annotation\n"))
#This is the header: "Gene name; A.; B.; E.; what_pattern; J.;J./n "

#Needs B. to E.
for mutation in mutations_dict:
    basename = mutation.split("_")
    if mutations_dict[mutation] in patterns:
        final_summary_file.write(str(goi_dict[basename[0]]+basename[1]+";"+pattern_dict[mutations_dict[mutation]] + ";Single SNP pattern;"+pfam_dict[gene]+"\n")) ## This wraps up the pattern recognition

#Needs F.
for gene in simplified_mutation_dict:
    if simplified_mutation_dict[gene] in patterns:
        final_summary_file.write(str(goi_dict[gene]+";"+ pattern_dict[simplified_mutation_dict[gene]] + ";all mutations added up;"+pfam_dict[gene]+"\n"))

#Needs G.
for gene in pres_abs_dict:
    if pres_abs_dict[gene] in patterns:
        final_summary_file.write(str(goi_dict[gene]+";"+ pattern_dict[pres_abs_dict[gene]] + ";Presence-Absence;"+pfam_dict[gene]+"\n")) ## This wraps up the pattern recognition

#Needs H.       
for gene in simple_pres_mutation_dict:
    if simple_pres_mutation_dict[gene] in patterns:
        final_summary_file.write(str(goi_dict[gene]+";"+ pattern_dict[simple_pres_mutation_dict[gene]] + ";Presence-Absence+simplified;"+pfam_dict[gene]+"\n"))

#Needs I.
for gene in pres_mutation_dict:
    if pres_mutation_dict[gene] in patterns:
        basename = gene.split(";")
        final_summary_file.write(str(goi_dict[basename[0]]+basename[1]+ ";"+ pattern_dict[pres_mutation_dict[gene]]+";Presence-Absence+detailed;"+pfam_dict[gene]+"\n"))

final_summary_file.close()