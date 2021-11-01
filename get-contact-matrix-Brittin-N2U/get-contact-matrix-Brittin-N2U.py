import mysql.connector as sql
from datetime import datetime
import getpass
import pandas as pd

##need to have mysql server running
## sudo systemctl start mysql
## see notes.txt on how to load databases from Brittin into server

##TODO get a robust version? probably can use imgNum as proxy for z-value

pw = getpass.getpass(prompt='Password: ', stream=None)
cnx = sql.connect(host='localhost', user='yinghong', password=pw, database='adult_db')
pw = 0

## adjacency2 column names:
## pre, post, preidx, postidx, imgNum, weight, preObj, postObj
query = 'SELECT pre, post, imgNum, weight FROM adjacency2;' ##add LIMIT 10 for test
cursor = cnx.cursor() #some iterator
cursor.execute(query)
contact = pd.DataFrame(cursor, columns = ['pre', 'post', 'imgNum', 'weight'])
cursor.close()
cnx.close()

#note: I checked that if A,B appears as pre,post, then B,A doesn't
#get total pixel of contact
contact_weight = contact[['pre','post','weight']]
contact_pixels = contact_weight.groupby(['pre','post']).sum().reset_index()
contact_pixels.columns = ['pre','post','pixels']
#get number of sections (if appear twice in an img, counts twice)
contact_imgNum = contact[['pre','post','imgNum']]
contact_imgNum = contact_imgNum.drop_duplicates()
contact_sections = contact_imgNum.groupby(['pre','post']).size().reset_index(name="sections")


###############################################
## save to file
#now = datetime.now()
#print("now =", now)
#dt_string = now.strftime("%d-%m-%Y-%H%M%S")

contact_pixels.to_csv('../data-aux/contact-Brittin-N2U-pixels.csv',encoding='utf-8-sig')
contact_sections.to_csv('../data-aux/contact-Brittin-N2U-sections.csv',encoding='utf-8-sig')

###############################################

#tested the list of cells found here against celllist from SI4
#some difference: 'HSNR/L', 'VC01' found in contact list
#but many more (354)found in celllist from SI4 than contact
