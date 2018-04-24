## import functiona
from Bio import Entrez, SeqIO
from ftplib import FTP
import ftplib
from pymongo import MongoClient
import glob
import os, sys, os.path
import re
import gzip
import time
## Today's date: September 19th
print time.strftime('%X %x %Z')

Entrez.tool = "September19"
Entrez.email = "snajar@uncc.edu"

ftpSite = "ftp.ncbi.nih.gov"
ftp = ftplib.FTP(ftpSite)
ftp.login()
ftp.cwd('/genbank/')

client = MongoClient()
mydb = client['test']
mycol = mydb['genbankDB']

dirFilenames = ftp.nlst()
# print dirFilenames
def sortingMethod(x):
	return(x[-11:])
sortedDirFilenames = sorted(dirFilenames, key = sortingMethod)  

def getAnnotation():
	with gzip.open(localPath+eachFile, mode='rb', compresslevel=9) as myWorkingFile:
		print "Reading: ", eachFile
		# get all sequence records for the specified genbank file
		records = [record for record in SeqIO.parse(myWorkingFile, "genbank")]
		print len(records)
		#print annotations for each sequence record
		for seqRecord in records:
			gid = seqRecord.id
			locus = seqRecord.name
			desc = seqRecord.description
			dblink = seqRecord.dbxrefs
			seq = seqRecord.seq
			# comments = seqRecord.annotations["comment"]
			source = seqRecord.annotations["source"]
			taxonomy = seqRecord.annotations["taxonomy"]
			pubDate = seqRecord.annotations["date"]
			organism = seqRecord.annotations["organism"]

			for ref in seqRecord.annotations["references"]:
				articleAuthors = ref.authors
				articleTitle = ref.title

			features = seqRecord.features
			country = "null"
			collection_date = "null"

			featSource = [featS for featS in features if featS.type == "source"]
			for featS in featSource:
				if 'country' in featS.qualifiers:
					country = featS.qualifiers['country']
				if 'collection_date' in featS.qualifiers:
					collection_date = featS.qualifiers['collection_date'][0]
					# collection_date = int(re.search(r'\d+', collection_date).group())
					# print "3.5\t" + str(collection_date)
				# else:
				# 	print "null"

			CDSproducts = []
			CDSgenes = []
			Gproducts = []
			Ggenes = []

			featsCDS = [featc for featc in features if featc.type == "CDS"]
			for featc in featsCDS:
				if 'product' in featc.qualifiers:
					productc=featc.qualifiers['product'][0].lower()
					CDSproducts.append(productc)
			 		# print featc
				elif 'gene' in featc.qualifiers:
					geneC = featc.qualifiers['gene'][0].lower()
					CDSgenes.append(geneC)

			featsGene = [featg for featg in features if featg.type == "gene"]
			for featg in featsGene:
				if 'product' in featg.qualifiers:
					product=featg.qualifiers['product'][0].lower()
					Gproducts.append(product)
				elif 'gene' in featg.qualifiers:
					gene = featg.qualifiers['gene'][0].lower()
					Ggenes.append(gene)

			# print gid + "\t" + source + "\t" + str(taxonomy) + "\t" + pubDate + "\t" + organism + "\t" + country + "\t" + str(CDSproducts) + "\t" + str(CDSgenes) + "\t" + str(Gproducts) + "\t" + str(Ggenes) + "\n"
			print "genbankID ************ " + gid

			mydb.test1.insert_one(
				{
				"id": gid,
				"type": "bacteria",
				"organism": organism,
				"locus": locus,
				"description": desc,
				"dblink": dblink,
				"taxonomy": taxonomy,
				"source": source,
				"country:": country,
				"publication date": pubDate,
				"collection_date": collection_date,
				"authors" : articleAuthors,
				"title": articleTitle,
				"sequence":	str(seq)},
				{"$set": {
					"CDSproducts": str(CDSproducts),
					"CDSgenes": str(CDSgenes),
					"Gproducts": str(Gproducts),
					"Ggenes": str(Ggenes),
				}}
				)
	
	# r'gbbct.*\.seq.gz$'

for eachFile in sortedDirFilenames:
	##################################### Pulling Bacteria Records ##########################################
	if re.match("gbbct245.seq.gz", eachFile):
		if not os.path.exists("/Users/snajar/Documents/pyworkspace/nuccoreFolder/bacteria/"):
			localPath = "/Users/snajar/Documents/pyworkspace/nuccoreFolder/bacteria/"
			os.makedirs(localPath, 0777)
			if os.path.isfile(localPath+eachFile):
				print "case 1:", eachFile, " exist."
				getAnnotation()
			else:
				print "file not exist, downloading...\ncase 2:", eachFile, "downloaded."
				ftp.retrbinary('RETR ' + eachFile, open(localPath+eachFile, 'wb').write)
				getAnnotation()
				 
		else:
			localPath = "/Users/snajar/Documents/pyworkspace/nuccoreFolder/bacteria/"
			if os.path.isfile(localPath+eachFile):
				print "case 3:", eachFile, "exist."
				getAnnotation()
			else:
				print "file not exist, downloading..."
				ftp.retrbinary('RETR ' + eachFile, open(localPath+eachFile, 'wb').write)
				print "case 4:", eachFile, "downloaded."
				getAnnotation()


	# print "1. Sleeping Mode... "
	# time.sleep(60)