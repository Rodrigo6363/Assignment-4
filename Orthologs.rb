###################################################
#                                                 #
#           Bioinformatic challenges              #
#             Rodrigo Álvarez Pardo               #
#                 Assignment 4                    #
#                                                 #
#                                                 #
###################################################

################# REQUIREMENTS ####################

require 'bio'
require 'stringio'
require 'io/console'

############### GLOBAL VARIABLES ##################

# References for the parameters used in BLAST :
#
# - https://www.ncbi.nlm.nih.gov/pubmed/18042555.
#   "Choosing BLAST options for better detection of orthologs as reciprocal best hits."
#   (Moreno-Hagelsieb G, Latimer K.)
#
# - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4094424/
#    "Quickly Finding Orthologs as Reciprocal Best Hits with BLAT, LAST, and UBLAST: How Much Do We Miss?"
#    (Natalie Ward and Gabriel Moreno-Hagelsieb)

$EVALUE = 10 ** -6
$OVERLAP = 50

################ CREATE OUTPUT FILE ####################

def output_file(file)    # Method for open the input file.
  
  if File.exists?(file)  # Delete in the case the file exits.
    File.delete(file)
  end
  
  return (File.open(file, "w"))
end

output = output_file ('output_ortologes.txt')

########## CREATING A DIRECTORY FOR DATABASES ##########

system ("mkdir Databases") # Build a directory that contains the databases

############### INPUT FILE AND FILE TYPE ###############

search_file = ("pep.fa") 
type_search_file = 'prot'
target_file = ("TAIR10_seq_20110103_representative_gene_model_updated")
type_target_file = 'nucl'

#################### DATABASE NAMES #####################

search_db = search_file.to_s + '_db'
target_db = target_file.to_s + '_db'


################### CREATING DATABASES ##################

system ( "makeblastdb -in #{search_file} -dbtype #{type_search_file} -out ./Databases/#{search_db}")
system (" makeblastdb -in #{target_file} -dbtype #{type_target_file} -out ./Databases/#{target_db}")

############### CREATING FACTORY OBJECTS ################
 
factory_search_file = Bio::Blast.local('blastx',"./Databases/#{search_db}")
# blastx uses a nucleotde sequence as input, and compares them whit a protein databases
factory_target_file = Bio::Blast.local('tblastn',"./Databases/#{target_db}")
# Tblastn compares a proteic sequence whit a nucleotide database
  
########## CREATE FASTA FORMAT FOR THE FILES #############

fasta_search_file = Bio::FastaFormat.open(search_file)
fasta_target_file = Bio::FastaFormat.open(target_file)

###################### REPORT BLAST #####################

target_hash = Hash.new                                                                              # We create a hash to store the target_file sequences
fasta_target_file.each do |target_sequence|
  target_hash[(target_sequence.entry_id).to_s] = (target_sequence.seq).to_s
end

count = 0                                                                                           # We establish a counter to know how many orthologs we have

fasta_search_file.each do |search_sequence|                                                         # Iterate over each sequence in the search file
  
  search_id = (search_sequence.entry_id).to_s                                                       # We store the id to know later if it is a reciprocal hit
  report_target = factory_target_file.query(search_sequence)
  
  
  if report_target.hits[0]
    
    target_id = (report_target.hits[0].definition.match(/(\w+\.\w+)|/)).to_s
    
    if (report_target.hits[0].evalue <= $EVALUE) and (report_target.hits[0].overlap >= $OVERLAP)    # We check parameters
      
      report_search = factory_search_file.query(">#{target_id}\n#{target_hash[target_id]}")
      # We look in the hash with the previous ID to get the sequence and query the factory
      
      if report_search.hits[0]
        
        match = (report_search.hits[0].definition.match(/(\w+\.\w+)|/)).to_s
        
        if (report_search.hits[0].evalue <= $EVALUE) and (report_search.hits[0].overlap >= $OVERLAP) # We check parameters
          
          if search_id == match                                                                      # If the match and search id match, this means that is a reciprocal hit
            
            output.puts "#{search_id}\t\t#{target_id}"                                               # We write the result in the output
            puts "#{search_id}\t#{target_id}"
            
            count += 1
          
          end
        end
      end
    end
  end
end

################### DELETE DATABASES ####################


system ("rm -r Databases") # We remove the directory at the end of the task

#########################################################

output.puts "\n\nNumber of orthologues found: #{count}"
puts "You can browse the output in the file output_ortologues.txt"

################# BONUS: NEXT STEPS #####################

# One of the things we can do to verify that our genes are
# orthologs is to conduct a phylogenetic study, building a
# tree in order to indetify them. This can be a bit tedious
# task because of the size of these files.

# We could also consult different online resources that are
# ortholog search engines or databases of these to contrast our results.




