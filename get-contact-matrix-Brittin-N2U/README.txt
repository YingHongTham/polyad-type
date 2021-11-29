######################################################################
get contact matrix from Brittin N2U
using sql

output to contact-Brittin-N2U-*.csv


for the .sql files, it stores information by issueing the commands
that create the database (CREATE DATABASE ...), table (CREATE TABLE ...), data (INSERT ...)
usually load to our sql server by
$ mysql -u yinghong -p < adult_databases.sql
but it seems the sql files are lacking the CREATE DATABASE commands...
manually added them myself:
(or just log in, create the empty database first, then do 
$ mysql -u yinghong -p adult_db < adult_databases.sql)

in adult_databases.sql:

CREATE DATABASE adult_db;
USE adult_db;

in l4_databases.sql:

CREATE DATABASE l4_db;
USE l4_db;

######################################################################
TODO: in hindsight, probably should have output in one table..

output two tables
contact-Brittin-N2U-sections.csv
contact-Brittin-N2U-pixels.csv

as edge list, i.e. has columns pre, post, sections/pixels

read with:
contact_sections = pd.read_csv('contact-Brittin-N2U-sections.csv', index_col=0)

######################################################################
