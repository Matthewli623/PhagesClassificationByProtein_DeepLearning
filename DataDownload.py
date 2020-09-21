from Bio import Entrez, SeqIO
from datetime import datetime
from multiprocessing import Process
from multiprocessing import Pool
import multiprocessing
import pandas, csv

# CSVs = ["Ino_Prot.csv","Myo_Prot.csv","Podo_Prot.csv","Sipho_Prot.csv"]
# SeqCSVs = ["Ino_Prot_Seq.csv","Myo_Prot_Seq.csv","Podo_Prot_Seq.csv","Sipho_Prot_Seq.csv"]
CSVs = ["Myo_Prot.csv", "Podo_Prot.csv", "Sipho_Prot.csv"]
SeqCSVs = ["Myo_Prot_Seq.csv", "Podo_Prot_Seq.csv", "Sipho_Prot_Seq.csv"]
dataL = []


def readCSV():
    for csv in CSVs:
        try:
            data = pandas.read_csv("/content/drive/My Drive/FYP/" + csv)
        except Exception as e:
            data = pandas.read_csv(csv)
        dataL.append(data)


def getProtSeqByID(id: str):
    count = 0
    while count < 4:
        try:
            print("start fetching {}".format(id))
            Entrez.email = "cyli7-c@my.cityu.edu.hk"
            handle = Entrez.esearch(db="protein", term=id)
            read_gene = Entrez.read(handle)
            # print(read_gene['IdList'])
            handle = Entrez.efetch(
                db="protein", id=read_gene["IdList"], rettype="fasta", retmode="text"
            )
            # print (handle.read())
            seq_record = SeqIO.read(handle, "fasta")
            # print (seq_record.seq)
            print("sucessfully fetched {}".format(id))
            return (seq_record.seq , seq_record.description)
        except Exception as e:
            # Try agin if any error
            print(e)
            count += 1
            print("{} is trying to fetch again".format(id))
    return ("","")


def writeSeq2CSV():
    for seqcsv, data in zip(SeqCSVs, dataL):
        with open(seqcsv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(
                ["prot_tag", "prot_id", "name", "size", "seq", "description"]
            )
            finishedCount = 0
            # Start loop
            PROCESS_NUM = 1
            while finishedCount != len(data):
                remainCount = len(data) - finishedCount
                processCount = 0
                # For create process
                if remainCount > PROCESS_NUM:
                    processCount = PROCESS_NUM
                else:
                    processCount = remainCount
                # each time processCount
                
                protTags = []
                protIds = []
                protNames = []
                protLengths = []
                processL = []
                
                for i in range(processCount):
                    protTags.append(data["prot_tag"][finishedCount + i])
                    protIds.append(data["prot_id"][finishedCount + i])
                    protNames.append(data["name"][finishedCount + i])
                    protLengths.append(data["size"][finishedCount + i])
                    now = datetime.now()
                    curTime = now.strftime("%H:%M:%S")
                    print(
                        "{} {} Doing index:{} {}".format(curTime, seqcsv, i, protIds[i])
                    )
                    finishedCount += 1
                pool = Pool(processCount)
                pool_outputs = pool.map(getProtSeqByID, protIds)

                for i in range(len(processL)):
                    seq = pool_outputs[i][0]
                    des = pool_outputs[i][1]
                    writer.writerow(
                        [
                            protTags[i],
                            protIds[i],
                            protNames[i],
                            protLengths[i],
                            seq,
                            des,
                        ]
                    )
            csvfile.truncate()
            csvfile.close()




if __name__=="__main__":
    

    readCSV()
    writeSeq2CSV()