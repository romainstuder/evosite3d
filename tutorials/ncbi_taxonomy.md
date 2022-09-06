Browsing the NCBI Taxonomy with Python

The Taxonomy browser is a fantastic tool to browser the Tree of Life:
http://www.ncbi.nlm.nih.gov/taxonomy
http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Root


I wrote a simple script in Python to browse it.

First, download the taxonomy archive in your folder:
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

Unpack it:
tar xvfz taxdump.tar.gz

It will output the following files:
citations.dmp
delnodes.dmp
division.dmp
gc.prt
gencode.dmp
merged.dmp
names.dmp
nodes.dmp
readme.txt


Only "names.dmp and "nodes.dmp" are important for us.

Download the script "parse_taxbrowser.py" from the GitHub and save it in the same folder:
https://github.com/romainstuder/evosite3d/blob/master/parse_taxbrowser.py


Then makes it executable:
chmod +x parse_taxbrowser.py

Then execute it in the same folder as "names.dmp and "nodes.dmp":
./parse_taxbrowser.py

And enjoy!

It will do different things:
1) Load the assocation between Tax_id and Names ("names.dmp")
2) Load the Tree of Life ("nodes.dmp")
3) Display the genealogy from root to Human
4) Find the common ancestor between Trichoplax and Human. 
5) Find descendants (internal or terminal leaves) of a specific taxa, based on division (i.e. "genus", "species").




