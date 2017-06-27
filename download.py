#Download genome data by species from PATRIC
#Clara Liu
#2017/06/19

from ftplib import FTP
import sys, os

#######if you have proxy

####fill in you proxy ip here
#site = FTP('1.1.1.1')

#site.set_debuglevel(1)
#msg = site.login('anonymous@ftp.patricbrc.org')

species = 'Mycobacterium_tuberculosis'
extension = '.ffn'

site = FTP("ftp.patricbrc.org")
site.login()
#site.cwd('/patric2/current_release/ec/')
site.cwd('/patric2/current_release/genomes_by_species/%s/' %species)


output = species+"/"+extension[1:]+"/"
if not os.path.exists(output):
    os.makedirs(output)

counter = 1
for filename in site.nlst():
  fhandle = open(output+filename+"."+extension, 'wb')
  print str(counter) + ': Getting ' + filename #for confort sake, shows the file that's being retrieved
  try:
    site.retrbinary('RETR ' + filename+'/'+filename+".PATRIC"+extension, fhandle.write)
    fhandle.close()
    counter += 1
  except:
    #os.remove(output+filename+"."+extension)
    print "%s not retrievable" %filename
