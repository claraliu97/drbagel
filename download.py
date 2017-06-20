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

site = FTP("ftp.patricbrc.org")
site.login()
#site.cwd('/patric2/current_release/ec/')
site.cwd('/patric2/current_release/genomes_by_species/%s/' %species)


output = species+"/"
os.mkdir(output)

counter = 1
for filename in site.nlst():
  fhandle = open(output+filename+'.fna', 'wb')
  print str(counter) + ': Getting ' + filename #for confort sake, shows the file that's being retrieved
  try:
    site.retrbinary('RETR ' + filename+'/'+filename+'.fna', fhandle.write)
    fhandle.close()
    counter += 1
  except:
    print "%s not retrievable" %filename