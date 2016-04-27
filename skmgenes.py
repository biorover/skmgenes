#!/usr/bin/env python
#Author: Sean Keane McKenzie, Graduate Fellow at The Rockefeller University, email: mckensk0@gmail.com
#requires custom library skmfunc (also by me)
import skmfunc
import StringIO
import copy
import string


def get_from_fasta(seq_name,fasta):
    if type(fasta) == file:
        fastaf=fasta
    else:
        try:
            fastaf=open(fasta)
        except:
            fastaf=StringIO.StringIO(fasta)
    seqline=fastaf.readline()
    locline=""
    while seqline != "":
        if ">" in seqline:
            bol=False
            if len(seq_name.split())== 1:
                if seqline[1:-1].split()[0] == seq_name:
                    bol=True
            else:
                if seqline[1:-1] == seq_name:
                    bol=True
            if bol:
                nextline=fastaf.readline()
                while nextline != "":
                    if ">" in nextline:
                        break
                    else:
                        locline=locline+nextline.split("\n")[0]
                        nextline=fastaf.readline()
                break
        seqline=fastaf.readline()
    fastaf.close()
    return locline

def grabseq(locus,cdslist,strand,genomelocation):
    """obtains gene sequence from a genome fasta file using scaffold/contg identifier and list of exon coordinates"""
    locline=get_from_fasta(locus,genomelocation)
    geneseq=""
    for i in cdslist:
        geneseq=geneseq+locline[i[0]-1:i[1]]
    if strand == "-":
        geneseq=skmfunc.revcomp(geneseq)
    return geneseq



class AnnotationDic(dict):
    """colection of genomic annotations"""
    #begin outdated
    # def remove_overlaps(self,metric="length"):
        # #build locus dictionary
        # locusdic=dict()
        # for annotation in self:
            # try:
                # locusdic[self[annotation].locus].append(annotation)
            # except:
                # locusdic[self[annotation].locus]=[annotation]
        # #checks for overlap
        # removelist=[]
        # for locus in locusdic:
            # #sort annotations by metric
            # metriclist=[]
            # for annotation in locusdic[locus]:
                # metriclist.append(eval("self[annotation]."+metric))
            # sort_loc=zip(*sorted(zip(metriclist,locusdic[locus])))[1]
            # for i in range(len(sort_loc)):
                # if not sort_loc[i] in removelist:
                    # for k in range(len(sort_loc))[i+1:]:
                        # if self[sort_loc[i]].stop < self[sort_loc[k]].start or self[sort_loc[i]].start > self[sort_loc[k]].stop:
                            # pass
                        # else:
                            # removelist.append(sort_loc[k])
        # for annotation in set(removelist):
            # del(self[annotation])
    # def write2gff(self,for_apollo=False,gff_format="Apollo",gff_object="gene"):
        # if for_apollo:
            # #seperate loci
            # loci_dic={}
            # for annotation in self:
                # try:
                    # loci_dic[self[annotation].locus].append(annotation)
                # except:
                    # loci_dic[self[annotation].locus]=[annotation]
        # else:
            # loci_dic={"pointer":self}
        # gfflist=[]
        # for locus in loci_dic:
            # locuslist=[]
            # for annotation in loci_dic[locus]:
                # locuslist.append(self[annotation].gff(gff_format=gff_format,gff_object=gff_object))
                # if locuslist[-1][-1] == "\n":
                    # locuslist[-1]=locuslist[-1][:-1]
            # if for_apollo:
                # locuslist.append("##FASTA\n>"+locus+"\n"+get_from_fasta(locus,self[annotation].genomelocation))
            # gfflist.append("\n".join(locuslist))
        # if for_apollo:
            # return gfflist
        # else:
            # return gfflist[0]
    #end outdated
    #adjusted for new annotation and exon classes
    def build_locus_dic(self):
        self.locus_dic={}
        for gene in self:
            locus=self[gene].locus
            if locus in self.locus_dic:
                self.locus_dic[locus].append(gene)
            else:
                self.locus_dic[locus]=[gene]
    def __add__(self,annother_annotation_dictionary):
        newdic=copy.deepcopy(self)
        for annotation in annother_annotation_dictionary:
            if annotation in self:
                print "annotation with ID "+annotation+" already in annotation dictionary"
            else:
                newdic[annotation]=copy.deepcopy(annother_annotation_dictionary[annotation])
        return newdic
    def rename(self,oldname,newname):
        self[oldname].rename(newname)
        self[newname]=self[oldname]
        del self[oldname]



class Exon(tuple):
    """exon object for annotations- really just a renamed tuple that can be assigned attributes"""
    pass


class Transcript(object):
    """A transcript object for populating the transcripts dictionary of an annotation object (skmgenes.Annotation)"""
    def __init__(self,ID="",name="",cdslist=[],exonlist=[]):
        self.cdslist=cdslist
        self.exonlist=exonlist
        self.name=name
        self.ID=ID
        self.update()

    def update(self):
        try:
            self.start=self.exonlist[0][0]
            self.stop=self.exonlist[-1][-1]
        except:
            self.start=0
            self.stop=0
        if self.name=="":
            self.name=self.ID

    def CDS2exons(self):
        self.exonlist=self.cdslist
        self.update()

    def exons2CDS(self):
        self.cdslist=self.exonlist
        self.update()

    def fillcdslist(self):
        for exon in self.exonlist:
            cdslist.append(exon.cds)

class Annotation(object):
    """genomic annotation"""
    def __init__(self,name="",ID="",locus="",strand="",transcripts={},genomelocation="null",source="PythonAnnotation"):
        self.name=name
        self.ID=ID
        self.locus=locus
        self.strand=strand
        self.genomelocation=genomelocation
        self.source=source
        self.transcripts={}
        self.update()
    
    def update(self):
        genestart=1000000000000000000
        genestop=0
        for transcript in self.transcripts:
            thandle=self.transcripts[transcript]
            thandle.update()
            if thandle.start < genestart:
                genestart=thandle.start
            if thandle.stop > genestop:
                genestop=thandle.stop
        self.start=genestart
        self.stop=genestop
        if self.name=="":
            self.name=self.ID
    
    def getseq(self):
        if self.genomelocation != "null":
            for transcript in transcripts:
                transcripts[transcript].CDS=grabseq(self.locus,transcripts[transcript].cdslist,self.strand,self.genomelocation)
                transcripts[transcript].cDNA=grabseq(self.locus,transcripts[transcript].exonlist,self.strand,self.genomelocation)
                seqlen=0
                for exon in transcripts[transcript].cdslist:
                    exonlen=exon[1]-exon[0]+1
                    seqlen=seqlen+exonlen
                transcripts[transcript].seqlength=seqlen

    
    def translate(self):
        """returns translation of annotation sequence"""
        if hasattr(self,"CDS"):
            return skmfunc.translatesequence(self.CDS)
        else:
            if self.genomelocation == "null":
                print "no genome location given"
            else:
                self.update()
                return skmfunc.translatesequence(self.CDS)
    
    def exonsequence(self,exon_number):
        cdslistcopy=self.cdslist[:]
        if self.strand == '-':
         cdslistcopy.reverse()
        priorlength=0
        for i in cdslistcopy[:exon_number-1]:
            priorlength=priorlength+i[1]-i[0]+1
        return self.sequence[priorlength:priorlength+cdslistcopy[exon_number-1][1]-cdslistcopy[exon_number-1][0]+1]
    
    def gff(self,gff_format="Apollo",gff_object="gene"):
        if gff_object == "gene":
            gfflines=["\t".join([self.locus,self.source,"mRNA",str(self.start),str(self.stop),".",self.strand,".","ID="+self.name+";\n"])]
            for cdsexon in self.cdslist:
                gfflines.append("\t".join([self.locus,self.source,"CDS",str(cdsexon[0]),str(cdsexon[1]),".",self.strand,".","Parent="+self.name+";\n"]))
            gff="".join(gfflines)
            if gff_format == "Apollo":
                gff=skmfunc.BGItoAPOLLO(gff)
        elif gff_object == "evidence":
            try:
                score=self.score
            except:
                score="."
            gfflines=["\t".join([self.locus,self.source,"match",str(self.start),str(self.stop),str(score),self.strand,".","ID="+self.name+";Name="+self.name+";\n"])]
            for cdsexon in self.cdslist:
                gfflines.append("\t".join([self.locus,self.source,"match_part",str(cdsexon[0]),str(cdsexon[1]),str(score),self.strand,".","ID="+self.name+":"+str(self.cdslist.index(cdsexon))+";Parent="+self.name+";\n"]))
            gff="".join(gfflines)
        return gff
    #fixed for new annotation and exon classes
    def clip(self,position,direction='+',level='nucleotide',remove_null_transcripts=True):
        null_transcripts=[]
        if not type(position) == tuple:
            if direction=='+':
                for transcript in self.transcripts:
                    removelist=[]
                    for i in range(len(self.transcripts[transcript].exonlist)):
                        exon=self.transcripts[transcript].exonlist[i]
                        if exon[1] >= position > exon[0]:
                            if level=='nucleotide':
                                newexondict=exon.__dict__.copy()
                                self.transcripts[transcript].exonlist[i]=Exon((exon[0],position-1))
                                self.transcripts[transcript].exonlist[i].__dict__=newexondict
                                if 'cds' in newexondict:
                                    if self.transcripts[transcript].exonlist[i].cds[1] >= position > self.transcripts[transcript].exonlist[i].cds[0]:
                                        self.transcripts[transcript].exonlist[i].cds=(self.transcripts[transcript].exonlist[i].cds[0],position-1)
                                    elif self.transcripts[transcript].exonlist[i].cds[0] >= position:
                                        del(self.transcripts[transcript].exonlist[i].cds)
                            elif level=='exon':
                                removelist.append(exon)
                        elif exon[0] >= position:
                                removelist.append(exon)
                    for exon in removelist:
                        self.transcripts[transcript].exonlist.remove(exon)
                    if len(self.transcripts[transcript].exonlist) == 0:
                        null_transcripts.append(transcript)
            elif direction=='-':
                for transcript in self.transcripts:
                    removelist=[]
                    for i in range(len(self.transcripts[transcript].exonlist)):
                        exon=self.transcripts[transcript].exonlist[i]
                        if exon[1] > position >= exon[0]:
                            if level=='nucleotide':
                                newexondict=exon.__dict__.copy()
                                self.transcripts[transcript].exonlist[i]=Exon((position+1,exon[1]))
                                self.transcripts[transcript].exonlist[i].__dict__=newexondict
                                if 'cds' in newexondict:
                                    if self.transcripts[transcript].exonlist[i].cds[1] > position >= self.transcripts[transcript].exonlist[i].cds[0]:
                                        self.transcripts[transcript].exonlist[i].cds=(position-1,self.transcripts[transcript].exonlist[i].cds[1])
                                    elif self.transcripts[transcript].exonlist[i].cds[1] <= position:
                                        del(self.transcripts[transcript].exonlist[i].cds)
                            elif level=='exon':
                                removelist.append(exon)
                        elif exon[1] <= position:
                            removelist.append(exon)
                    for exon in removelist:
                        self.transcripts[transcript].exonlist.remove(exon)
                    if len(self.transcripts[transcript].exonlist) == 0:
                        null_transcripts.append(transcript)
            if remove_null_transcripts:
                for transcript in null_transcripts:
                    del(self.transcripts[transcript])
        else:
            for transcript in self.transcripts:
                removelist=[]
                for i in range(len(self.transcripts[transcript].exonlist)):
                    exon=self.transcripts[transcript].exonlist[i]
                    if exon[0] <= position[0] <= exon[1] <= position[1]:
                        if level=='nucleotide':
                            newexondict=exon.__dict__.copy()
                            self.transcripts[transcript].exonlist[i]=Exon((exon[0],position[0]-1))
                            self.transcripts[transcript].exonlist[i].__dict__=newexondict
                            if 'cds' in newexondict:
                                if self.transcripts[transcript].exonlist[i].cds[1] >= position[0] > self.transcripts[transcript].exonlist[i].cds[0]:
                                    self.transcripts[transcript].exonlist[i].cds=(self.transcripts[transcript].exonlist[i].cds[0],position[0]-1)
                                elif self.transcripts[transcript].exonlist[i].cds[0] >= position[0]:
                                    del(self.transcripts[transcript].exonlist[i].cds)
                        elif level=='exon':
                            removelist.append(exon)
                    elif position[0] <= exon[0] <= position[1] <= exon[1]:
                        if level=='nucleotide':
                            newexondict=exon.__dict__.copy()
                            self.transcripts[transcript].exonlist[i]=Exon((position[1]+1,exon[1]))
                            self.transcripts[transcript].exonlist[i].__dict__=newexondict
                            if 'cds' in newexondict:
                                if self.transcripts[transcript].exonlist[i].cds[1] >= position[1] > self.transcripts[transcript].exonlist[i].cds[0]:
                                    self.transcripts[transcript].exonlist[i].cds=(position[1]+1,self.transcripts[transcript].exonlist[i].cds[1])
                                elif self.transcripts[transcript].exonlist[i].cds[1] <= position[1]:
                                    del(self.transcripts[transcript].exonlist[i].cds)
                        elif level=='exon':
                            removelist.append(exon)
                    elif position[0] >= exon[0] >= exon[1] >= position[1]:
                            removelist.append(exon)
                for exon in removelist:
                    self.transcripts[transcript].exonlist.remove(exon)
                if len(self.transcripts[transcript].exonlist) == 0:
                    null_transcripts.append(transcript)
        self.update()

    def rename(self,newname,rename_transcripts=True):
        self.name=newname
        self.ID=newname
        if rename_transcripts==True:
            count=0
            deletelist=[]
            for transcript in copy.copy(self.transcripts):
                self.transcripts[newname+'-R'+string.uppercase[count]]=self.transcripts[transcript]
                deletelist.append(transcript)
                self.transcripts[newname+'-R'+string.uppercase[count]].name=newname+'-R'+string.uppercase[count]
                self.transcripts[newname+'-R'+string.uppercase[count]].ID=newname+'-R'+string.uppercase[count]
            for i in deletelist:
                del self.transcripts[i]

##begin outdated



def oldgff2annotations(gff,genomelocation=False):
    """reads annotations from a gff file or string and returns a dictionary of python annotation objects"""
    if type(gff) == "file":
        gfflines=gff.read().split(">")[0].split("\n")
        if len(gff.read().split(">")) > 1:
            gfflocation=gff.name
        else:
            gfflocation="null"
    else:
        try:
            gfflines=open(gff).read().split(">")[0].split("\n")
            if len(open(gff).read().split(">")) > 1:
                gfflocation=gff
            else:
                gfflocation="null"
        except:
            gfflines=gff.split(">")[0].split("\n")
            gfflocation="null"
    if genomelocation:
        gfflocation=genomelocation
    annotation_dic=AnnotationDic()
    for line in gfflines:
        if line == "":
            pass
        elif line[0] == "#":
            pass
        elif line.split("\t")[2].lower() == "mrna":
            cdslines=[]
            exonlines=[]
            fields=line.split("\t")
            genename=fields[-1].split(";")[0].split("=")[1]
            strand=fields[6]
            locus=fields[0]
            line_index=gfflines.index(line)+1
            while line_index < len(gfflines):
                fields=gfflines[line_index].split("\t")
                if fields == [""]:
                    pass
                elif fields[0][0] == "#":
                    pass
                elif fields[2].lower() == "cds":
                    cdsindex1=int(fields[3])
                    cdsindex2=int(fields[4])
                    if cdsindex1 < cdsindex2:
                        cdslines.append((cdsindex1,cdsindex2))
                    else:
                        cdslines.append((cdsindex2,cdsindex1))
                elif fields[2].lower() == "exon":
                    exonindex1=int(fields[3])
                    exonindex2=int(fields[4])
                    if exonindex1 < exonindex2:
                        exonlines.append((exonindex1,exonindex2))
                    else:
                        exonlines.append((exonindex2,exonindex1))
                elif fields[2].lower() == "mrna":
                    break
                line_index = line_index + 1
            if cdslines != []:
                if cdslines[0] > cdslines[-1]:
                    cdslines.reverse()
            if exonlines != []:
                if exonlines[0] > exonlines[-1]:
                    exonlines.reverse()
            annotation_dic[genename]=Annotation(name=genename,locus=locus,cdslist=cdslines,strand=strand,exonlist=exonlines,genomelocation=gfflocation)
    return annotation_dic



###end outdated







#fixed for new annotation class and exon class


def gff3_2_annotations(gff,genomelocation="null",fillCDS=False,trimUTR=False):
    """reads annotations from a gff file or string and returns a dictionary of python annotation objects"""
    gfflines=gff.replace('\r','').split('\n')
    annotation_dic=AnnotationDic()
    transcripttracker={}
    while "" in gfflines:
        gfflines.remove("")
    for line in gfflines:
        if line[0]!='#' and len(line.split('\t')) > 6:
            features=line.split('\t')
            if features[2]=='gene':
                deflinedic={}
                defline=features[8].replace('Name','name').split(';')
                for i in defline:
                    deflinedic[i.split('=')[0]]=i.split('=')[1]
                geneid=deflinedic['ID']
                if 'name' in deflinedic:
                    genename=deflinedic['name']
                else: genename=geneid[:]
                annotation_dic[geneid]=Annotation(ID=geneid,name=genename,locus=features[0],strand=features[6],source=features[1],genomelocation=genomelocation)
                for i in deflinedic:
                    if i != 'name' and i != 'ID':
                        setattr(annotation_dic[geneid],i,deflinedic[i])
            elif features[2]=='mRNA' or features[2]=='transcript' or features[2]=='match':
                deflinedic={}
                defline=features[8].replace('Name','name').split(';')
                for i in defline:
                    if i != "":
                        deflinedic[i.split('=')[0]]=i.split('=')[1]
                if 'Parent' in deflinedic: geneid=deflinedic['Parent']
                else: geneid=deflinedic['ID']+'-gene'
                transcriptid=deflinedic['ID']
                if 'name' in deflinedic: transcriptname=deflinedic['name']
                else: transcriptname=transcriptid[:]
                if geneid in annotation_dic:
                    annotation_dic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=[],cdslist=[])
                    transcripttracker[transcriptid]=geneid
                else:
                    annotation_dic[geneid]=Annotation(ID=geneid,name=transcriptname,locus=features[0],strand=features[6],source=features[1],genomelocation=genomelocation)
                    annotation_dic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=[],cdslist=[])
                    annotation_dic[geneid].update()
                    transcripttracker[transcriptid]=geneid
                if features[2]=='match':
                    annotation_dic[geneid].type='Evidence'
                for i in deflinedic:
                    if i != 'name' and i != 'ID':
                        setattr(annotation_dic[geneid].transcripts[transcriptid],i,deflinedic[i])
            elif features[2]=='exon' and not trimUTR or features[2]=='CDS' and trimUTR or features[2]=='match_part':
                deflinedic={}
                defline=features[8].replace('Name','name').split(';')
                for i in defline:
                    deflinedic[i.split('=')[0]]=i.split('=')[1]
                transcriptid=deflinedic['Parent']
                exoncoords=(int(features[3]),int(features[4]))
                if exoncoords[0] > exoncoords[1]:
                    exoncoords=list(exoncoords)
                    exoncoords.reverse()
                    exoncoords=tuple(exoncoords)
                if transcriptid in transcripttracker:
                    geneid=transcripttracker[transcriptid]
                elif transcriptid in annotation_dic:
                    geneid=transcriptid
                    annotation_dic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=[],cdslist=[])
                    transcripttracker[transcriptid]=geneid
                else:
                    geneid=transcriptid
                    annotation_dic[geneid]=Annotation(ID=geneid,name=geneid,locus=features[0],strand=features[6],source=features[1],genomelocation=genomelocation)
                    annotation_dic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=[],cdslist=[])
                    transcripttracker[transcriptid]=geneid
                annotation_dic[geneid].transcripts[transcriptid].exonlist.append(Exon(exoncoords))
                if fillCDS:
                    annotation_dic[geneid].transcripts[transcriptid].exonlist[-1].cds=exoncoords
                    #remove when CDS list depreciated
                    annotation_dic[geneid].transcripts[transcriptid].cdslist.append(exoncoords)
                    #
                try: annotation_dic[geneid].transcripts[transcriptid].exonlist[-1].score=float(features[5])
                except: pass
                annotation_dic[geneid].transcripts[transcriptid].exonlist.sort()
                annotation_dic[geneid].transcripts[transcriptid].update()
                annotation_dic[geneid].update()
            elif features[2]=='CDS' and not trimUTR and not fillCDS:
                deflinedic={}
                defline=features[8].replace('Name','name').split(';')
                for i in defline:
                    if i != "":
                        deflinedic[i.split('=')[0]]=i.split('=')[1]
                transcriptid=deflinedic['Parent']
                exoncoords=(int(features[3]),int(features[4]))
                if exoncoords[0] > exoncoords[1]:
                    exoncoords=list(exoncoords)
                    exoncoords.reverse()
                    exoncoords=tuple(exoncoords)
                if transcriptid in transcripttracker:
                    geneid=transcripttracker[transcriptid]
                elif transcriptid in annotation_dic:
                    geneid=transcriptid
                    annotation_dic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=[],cdslist=[])
                    transcripttracker[transcriptid]=geneid
                else:
                    geneid=transcriptid
                    annotation_dic[geneid]=Annotation(ID=geneid,name=geneid,locus=features[0],strand=features[6],source=features[1],genomelocation=genomelocation)
                    annotation_dic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid)
                    transcripttracker[transcriptid]=geneid
                foundexon=False
                partialoverlap=False
                for exon in annotation_dic[geneid].transcripts[transcriptid].exonlist:
                    if exon[0] <= exoncoords[0] <= exoncoords[1] <= exon[1]:
                        exon.cds=exoncoords
                        foundexon=True
                        partialoverlap=False
                        break
                    elif exoncoords[0] > exon[1]:
                        pass
                    elif exoncoords[1] < exon[0]:
                        pass
                    else:
                        partialoverlap=True
                        foundexon=True
                        exon.cds=exoncoords
                if partialoverlap:
                    print 'error, gene '+geneid+' transcript '+transcriptid+' CDS '+str(exoncoords)+' not completely contained within exon'+str(exon)
                if not foundexon:
                    annotation_dic[geneid].transcripts[transcriptid].exonlist.append(Exon(exoncoords))
                    annotation_dic[geneid].transcripts[transcriptid].exonlist[-1].cds=exoncoords
                    annotation_dic[geneid].transcripts[transcriptid].exonlist.sort()
                    annotation_dic[geneid].transcripts[transcriptid].update()
                    annotation_dic[geneid].update()
                #remove when CDS list depreciated
                annotation_dic[geneid].transcripts[transcriptid].cdslist.append(exoncoords)
                #
    return annotation_dic


def cegmagff2annotations(cegmagff,genomelocation="null",fillCDS=True):
    """converts gff output from CEGMA (core eukaryotic gene mapping analysis) to skmgenes annotation format"""
    annotation_dic=AnnotationDic()
    gfflines=cegmagff.split('\n')
    for line in gfflines:
        features=line.split('\t')
        if len(features) > 6:
            geneid=features[8]
            if geneid in annotation_dic:
                pass
            else:
                annotation_dic[geneid]=Annotation(ID=geneid,name=geneid,locus=features[0],strand=features[6],source=features[1],genomelocation=genomelocation)
                annotation_dic[geneid].transcripts[geneid+'-RA']=Transcript(ID=geneid+'-RA',exonlist=[],cdslist=[])
            annotation_dic[geneid].transcripts[geneid+'-RA'].exonlist.append(Exon((int(features[3]),int(features[4]))))
            if fillCDS:
                annotation_dic[geneid].transcripts[geneid+'-RA'].cdslist.append(Exon((int(features[3]),int(features[4]))))
                annotation_dic[geneid].transcripts[geneid+'-RA'].exonlist[-1].cds=(int(features[3]),int(features[4]))
                annotation_dic[geneid].transcripts[geneid+'-RA'].cdslist.sort()
            try: annotation_dic[geneid].transcripts[geneid+'-RA'].exonlist[-1].score=float(features[5])
            except: pass
            annotation_dic[geneid].transcripts[geneid+'-RA'].exonlist.sort()
            annotation_dic[geneid].transcripts[geneid+'-RA'].update()
            annotation_dic[geneid].update()
    return annotation_dic


def set_difference(annotation_dic_to_clip,annotation_dic_cookie_cutter,level='exon',reports=False, stranded=True, ignoreexonless=True):
    """returns a copy of an annotation dictionary (annotation_dic_to_clip) with any overlap with a second annotation dictionary removed"""
    newandic=copy.deepcopy(annotation_dic_to_clip)
    newandic.build_locus_dic()
    cc=[]
    if level == 'gene':
        for i in annotation_dic_cookie_cutter:
            cc.append((annotation_dic_cookie_cutter[i].start,annotation_dic_cookie_cutter[i].stop,annotation_dic_cookie_cutter[i].locus,annotation_dic_cookie_cutter[i].strand,i))
    elif level == 'exon' or level == 'nucleotide':
        for i in annotation_dic_cookie_cutter:
            for transcript in annotation_dic_cookie_cutter[i].transcripts:
                for exon in annotation_dic_cookie_cutter[i].transcripts[transcript].exonlist:
                    cc.append((exon[0],exon[1],annotation_dic_cookie_cutter[i].locus,annotation_dic_cookie_cutter[i].strand,i))
    newgenecount=1
    removelist=[]
    cliplist=[]
    clippedbydic={}
    addlist=[]
    for cutline in cc:
        cutcoords=(cutline[0],cutline[1])
        locus=cutline[2]
        strand=cutline[3]
        cutname=cutline[4]
        if locus in newandic.locus_dic:
            for tgeneid in newandic.locus_dic[locus]:
                tgene=newandic[tgeneid]
                evalexception=False
                if stranded:
                    if tgene.strand != strand:
                        evalexception=True
                if ignoreexonless:
                    if tgene.stop==0:
                        evalexception=True
                if tgene.start <= cutcoords[0] <= tgene.stop <= cutcoords[1] and not evalexception:
                    if level == 'exon' or level == 'nucleotide':
                        tgene.clip(cutcoords[0],direction='+',level=level)
                        cliplist.append(tgeneid)
                        clippedbydic[tgeneid]=cutname
                    elif level == 'gene':
                        removelist.append(tgeneid)
                elif cutcoords[0] <= tgene.start <= tgene.stop <= cutcoords[1] and not evalexception:
                    removelist.append(tgeneid)
                elif cutcoords[0] <= tgene.start <= cutcoords[1] <= tgene.stop and not evalexception:
                    if level == 'exon' or level == 'nucleotide':
                        tgene.clip(cutcoords[1],direction='-',level=level)
                        cliplist.append(tgeneid)
                        clippedbydic[tgeneid]=cutname
                    elif level == 'gene':
                        removelist.append(tgeneid)
                elif tgene.start <= cutcoords[0] <= cutcoords[1] <= tgene.stop and not evalexception:
                    if level == 'exon' or level == 'nucleotide':
                        tgene.clip((cutcoords[0],cutcoords[1]),direction='+', level=level)
                        cliplist.append(tgeneid)
                        clippedbydic[tgeneid]=cutname
                    elif level == 'gene':
                        removelist.append(tgeneid)
                tgene.update()
                if tgene.stop==0 and not evalexception:
                    removelist.append(tgeneid)
                    try: cliplist.remove(tgeneid)
                    except: pass
    cliplist=list(set(cliplist))
    for tgeneid in set(removelist):
        del newandic[tgeneid]
        try:
            cliplist.remove(tgeneid)
        except: pass
    if reports:
        return (newandic,list(set(removelist)),cliplist,clippedbydic)
    else:
        return newandic

def exonerate2annotations(exonerate_file,genomelocation="null",fillCDS=False):
    if type(exonerate_file) == file:
        exoneratelist=exonerate_file.read().split("C4 Alignment:")
    else:
        try:
            exoneratelist=open(exonerate_file).read().split("C4 Alignment:")
        except:
            exoneratelist=exonerate_file.split("C4 Alignment:")
    annotation_dic=AnnotationDic()
    for alignment in exoneratelist:
        for line in alignment.split("\n"):
            if line[:7] == "vulgar:":
                vulgar=line.split()
                genename=vulgar[1]+"_"+vulgar[5]+":"+vulgar[6]+"-"+vulgar[7]
                locus=vulgar[5]
                intron=0
                strand=vulgar[8]
                cdslistlist=[[int(vulgar[6]),eval(vulgar[6]+"-"+strand+"1")]]
                if strand == "+":
                    cdslistlist[0][0]=cdslistlist[0][0]+1
                    cdslistlist[0][1]=cdslistlist[0][1]+1
                for i in range(len(vulgar))[10:]:
                    if i%3 == 1:
                        if vulgar[i] == "M" or vulgar[i] == "G" or vulgar[i] == "S" or vulgar[i] == "F":
                            cdslistlist[-1][1]=eval(str(cdslistlist[-1][1])+strand+vulgar[i+2])
                        elif vulgar[i] == "5":
                            cdslistlist.append([eval(str(cdslistlist[-1][1])+strand+"2"),0])
                        elif vulgar[i] == "I":
                            cdslistlist[-1][0]=eval(str(cdslistlist[-1][0])+strand+vulgar[i+2])
                        elif vulgar[i] == "3":
                            cdslistlist[-1][0]=eval(str(cdslistlist[-1][0])+strand+str(3))
                            cdslistlist[-1][1]=eval(str(cdslistlist[-1][0])+"-"+strand+str(1))
                exonlist=[]
                for i in cdslistlist:
                    if strand == "+":
                        exonlist.append(Exon((i[0],i[1])))
                    else:
                        exonlist.append(Exon((i[1],i[0])))
                if strand == "-":
                    exonlist.reverse()
                cdslist=[]
                if fillCDS:
                    for exon in exonlist:
                        exon.cds=tuple(exon)
                        cdslist.append(tuple(exon))
                annotation_dic[genename]=Annotation(ID=genename,name=genename,locus=locus,strand=strand,source="Exonerate",genomelocation=genomelocation)
                annotation_dic[genename].type='evidence'
                annotation_dic[genename].score=vulgar[9]
                annotation_dic[genename].qstart=vulgar[2]
                annotation_dic[genename].qend=vulgar[3]
                annotation_dic[genename].transcripts[genename+'-RA']=Transcript(ID=genename+'-RA',name=genename+'-RA',exonlist=exonlist,cdslist=cdslist)
                annotation_dic[genename].update()
    return annotation_dic





#fixed for new annotation class


def PASAgtf2annotations(PASAgtf,genomelocation="null"):
    """converts GTF output from PASA (program to assemble spliced allignments) to skmgenes annotation dictionary format"""
    gtflines=PASAgtf.split('\n')
    PASAdic=AnnotationDic()
    for lineindex in range(len(gtflines)):
        line=gtflines[lineindex]
        features=line.split('\t')
        if len(features) > 6:
            if features[2]=='transcript':
                indi=lineindex+1
                source=features[1]
                locus=features[0]
                strand=features[6]
                geneid=features[8].split('"')[1]
                transcriptid=features[8].split('"')[3]
                exonlist=[]
                while indi < len(gtflines) and gtflines[indi].split('\t')[2] != 'transcript':
                    newfeatures=gtflines[indi].split('\t')
                    if len(newfeatures) > 6:
                        if newfeatures[2]=='exon':
                            exonlist.append((int(newfeatures[3]),int(newfeatures[4])))
                    indi=indi+1
                if not geneid in PASAdic:
                    PASAdic[geneid]=Annotation(ID=geneid,name=geneid,locus=locus,strand=strand,source=source,genomelocation=genomelocation)
                PASAdic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=exonlist)
                PASAdic[geneid].update()
    return PASAdic

def MAKERevidence_gff2annotations(MAKERevidence_gff,genomelocation="null",populateCDS=False):
    """converts GFF output from MAKER to skmgenes annotation dictionary format"""
    gfflines=MAKERevidence_gff.split('\n')
    MAKERdic=AnnotationDic()
    while '' in gfflines:
        gfflines.remove('')
    for lineindex in range(len(gfflines)):
        line=gfflines[lineindex]
        features=line.split('\t')
        if len(features) > 6:
            if features[2] != 'match_part':
                indi=lineindex+1
                source=features[1]
                locus=features[0]
                strand=features[6]
                geneid=features[8].split(';')[0][3:]
                name=features[8].split(';')[1][5:]
                try: score=float(features[5])
                except: pass
                target=name[:]
                evtype=features[2]
                transcriptid=geneid+'onlyTranscript'
                exonlist=[]
                while indi < len(gfflines) and gfflines[indi].split('\t')[2]=='match_part':
                    newfeatures=gfflines[indi].split('\t')
                    if len(newfeatures) > 6:
                        exonlist.append((int(newfeatures[3]),int(newfeatures[4])))
                    indi=indi+1
                if len(exonlist) > 1:
                    if exonlist[0][0] > exonlist[1][0]:
                        exonlist.reverse()
                MAKERdic[geneid]=Annotation(ID=geneid,name=name,locus=locus,strand=strand,source=source,genomelocation=genomelocation)
                try: MAKERdic[geneid].score=score
                except: pass
                MAKERdic[geneid].target=target
                MAKERdic[geneid].evtype=evtype
                MAKERdic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=exonlist)
                if populateCDS:
                    MAKERdic[geneid].transcripts[transcriptid].exons2CDS()
                MAKERdic[geneid].update()
    return MAKERdic

def augustus_gff2annotations(augustus_gff,genomelocation="null"):
    gfflines=augustus_gff.split('\n')
    AUGdic=AnnotationDic()
    for lineindex in range(len(gfflines)):
        line=gfflines[lineindex]
        features=line.split('\t')
        if len(features) > 6:
            if features[2]=='transcript':
                indi=lineindex+1
                source=features[1]
                locus=features[0]
                strand=features[6]
                geneid=features[8].split('.')[0]
                name=geneid[:]
                transcriptid=features[8]
                score=float(features[5])
                exonlist=[]
                while gfflines[indi][0]!='#' and len(gfflines[indi].split('\t')) > 6:
                    newfeatures=gfflines[indi].split('\t')
                    if newfeatures[2]=='CDS':
                        exonlist.append((int(newfeatures[3]),int(newfeatures[4])))
                    indi=indi+1
                AUGdic[geneid]=Annotation(ID=geneid,name=name,locus=locus,strand=strand,source=source,genomelocation=genomelocation)
                AUGdic[geneid].score=score
                AUGdic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=exonlist,cdslist=exonlist)
                AUGdic[geneid].update()
    return AUGdic


def genemark_gtf2annotations(genemark_gtf,genomelocation="null"):
    gtflines=genemark_gtf.split('\n')
    GMdic=AnnotationDic()
    for line in gtflines:
        features=line.split('\t')
        if len(features) > 6:
            if features[2]=='CDS':
                geneid=features[8].split('"')[1]
                transcriptid=geneid+'_onlyTranscript'
                if geneid in GMdic:
                    GMdic[geneid].transcripts[transcriptid].exonlist.append((int(features[3]),int(features[4])))
                    GMdic[geneid].transcripts[transcriptid].cdslist.append((int(features[3]),int(features[4])))
                    GMdic[geneid].transcripts[transcriptid].update()
                    GMdic[geneid].update()
                else:
                    GMdic[geneid]=Annotation(ID=geneid,name=geneid,locus=features[0],strand=features[6],source=features[1],genomelocation=genomelocation)
                    GMdic[geneid].transcripts[transcriptid]=Transcript(ID=transcriptid,exonlist=[(int(features[3]),int(features[4]))],cdslist=[(int(features[3]),int(features[4]))])
                    GMdic[geneid].update()
    return GMdic


def write2EVMprediction_gff(annotation_dic):
    outlines=[]
    for gene in annotation_dic:
        ghandle=annotation_dic[gene]
        outlines.append('\t'.join([ghandle.locus,ghandle.source,'gene',str(ghandle.start),str(ghandle.stop),'.',ghandle.strand,'.','ID='+ghandle.ID+';Name='+ghandle.name]))
        for transcript in ghandle.transcripts:
            thandle=ghandle.transcripts[transcript]
            outlines.append('\t'.join([ghandle.locus,ghandle.source,'mRNA',str(thandle.start),str(thandle.stop),'.',ghandle.strand,'.','ID='+thandle.ID+';Parent='+ghandle.ID]))
            for i in range(len(thandle.cdslist)):
                outlines.append('\t'.join([ghandle.locus,ghandle.source,'exon',str(thandle.cdslist[i][0]),str(thandle.cdslist[i][1]),'.',ghandle.strand,'.','ID='+thandle.ID+'_exon'+str(i)+';Parent='+thandle.ID]))
                outlines.append('\t'.join([ghandle.locus,ghandle.source,'CDS',str(thandle.cdslist[i][0]),str(thandle.cdslist[i][1]),'.',ghandle.strand,'.','ID='+thandle.ID+'_CDS'+str(i)+';Parent='+thandle.ID]))
    return '\n'.join(outlines)

def write2EVMevidence_gff(annotation_dic,evtype='cDNA_match'):
    outlines=[]
    for gene in annotation_dic:
        ghandle=annotation_dic[gene]
        for transcript in ghandle.transcripts:
            thandle=ghandle.transcripts[transcript]
            if hasattr(thandle,'score'):
                tscore=str(thandle.score)
            else:
                tscore='.'
            if hasattr(ghandle,'target'):
                ttarget=ghandle.target
            elif hasattr(thandle,'target'):
                ttarget=thandle.target
            else:
                ttarget=ghandle.name
            if hasattr(ghandle,'evtype'):
                tevtype=ghandle.evtype
            else:
                tevtype=evtype
            for exon in thandle.exonlist:
                outlines.append('\t'.join([ghandle.locus,ghandle.source,tevtype,str(exon[0]),str(exon[1]),tscore,ghandle.strand,'.','ID='+ghandle.ID+';Target='+ttarget]))
    return '\n'.join(outlines)

#adjusted for new annotation and exon classes

def write2fasta(annotation_dic,genomelocation='defualt',seqtype="nucl",feature="CDS",findframe=False,split_name=False):
    """returns a fasta format string containing sequences for all transcripts in all annotations in an annotation dictionary, seqtype can be "nucl" or "prot", feature can be "CDS" or "exon", findframe=True executes skmfunc.findframe on each transcript sequence"""
    if genomelocation=='default':
        indi=list(annotation_dic)[0]
        try:
            genomefile=open(annotation_dic[indi].genomelocation)
        except:
            print "Cannot open genome file "+annotation_dic[indi].genomelocation
    else:
        genomefile=open(genomelocation)
    genome=genomefile.read().split('>')[1:]
    genomedic={}
    for i in genome:
        name=i.split('\n')[0]
        if split_name:
            name=name.split()[0]
        seq=''.join(i.split('\n')[1:])
        genomedic[name]=seq
    fastalist=[]
    for gene in annotation_dic:
        for transcript in annotation_dic[gene].transcripts:
            transcriptseq=""
            for exon in annotation_dic[gene].transcripts[transcript].exonlist:
                locline=genomedic[annotation_dic[gene].locus]
                transcriptseq=transcriptseq+locline[exon.cds[0]-1:exon.cds[1]]
            if annotation_dic[gene].strand == "-":
                transcriptseq=skmfunc.revcomp(transcriptseq)
            if findframe:
                transcriptseq=skmfunc.findframe(transcriptseq)
            if seqtype=='prot':
                transcriptseq=skmfunc.translatesequence(transcriptseq)
            fastalist.append('>'+transcript+'\n'+transcriptseq)
    fasta='\n'.join(fastalist)
    return fasta

def annotation2gff3(annotation,is_evidence='auto',supress_defline_info=False):
    annotation.update()
    outlines=[]
    isev=False
    if 'type' in annotation.__dict__:
        if annotation.type=='evidence':
            isev=True
    if is_evidence=='yes':
        isev=True
    elif is_evidence=='no':
        isev=False
    try: gscore=str(annotation.score)
    except: gscore='.'
    deflist=['ID='+annotation.ID,'Name='+annotation.name]
    for i in annotation.__dict__:
        if i != 'score' and i != 'name' and i != 'ID' and i != 'strand' and i != 'start' and i != 'stop' and i != 'locus' and i != 'genomelocation' and i != 'transcripts' and i != 'source' and not supress_defline_info:
            deflist.append(i+'='+str(annotation.__dict__[i]))
    if not isev:
        outlines.append('\t'.join([annotation.locus,annotation.source,'gene',str(annotation.start),str(annotation.stop),gscore,annotation.strand,'.',';'.join(deflist)]))
    for transcript in annotation.transcripts:
        thandle=annotation.transcripts[transcript]
        exonlines=[]
        cdslines=[]
        deflist=['ID='+thandle.ID,'Name='+thandle.name,'Parent='+annotation.ID]
        try: tscore=str(thandle.score)
        except:
            tscore=gscore
        for i in annotation.__dict__:
            if i != 'score' and i != 'name' and i != 'ID' and i != 'strand' and i != 'start' and i != 'stop' and i != 'locus' and i != 'genomelocation' and i != 'transcripts' and i != 'source' and not supress_defline_info:
                deflist.append(i+'='+str(annotation.__dict__[i]))
        if isev:
            outlines.append('\t'.join([annotation.locus,annotation.source,'match',str(thandle.start),str(thandle.stop),tscore,annotation.strand,'.',';'.join(deflist)]))
        else:
            outlines.append('\t'.join([annotation.locus,annotation.source,'mRNA',str(thandle.start),str(thandle.stop),tscore,annotation.strand,'.',';'.join(deflist)]))
        for i in range(len(thandle.exonlist)):
            exon=thandle.exonlist[i]
            try: escore=str(exon.score)
            except: escore=tscore
            if isev:
                exonlines.append('\t'.join([annotation.locus,annotation.source,'match_part',str(exon[0]),str(exon[1]),escore,annotation.strand,'.','ID='+thandle.ID+':match_part'+str(i)+';Parent='+thandle.ID]))
            else:
                exonlines.append('\t'.join([annotation.locus,annotation.source,'exon',str(exon[0]),str(exon[1]),escore,annotation.strand,'.','ID='+thandle.ID+':exon'+str(i)+';Parent='+thandle.ID]))
                if 'cds' in exon.__dict__:
                    cdslines.append('\t'.join([annotation.locus,annotation.source,'CDS',str(exon[0]),str(exon[1]),escore,annotation.strand,'.','ID='+thandle.ID+':exon'+str(i)+'-CDS;Parent='+thandle.ID]))
        outlines=outlines+exonlines+cdslines
    return '\n'.join(outlines)

def write2cufflinks_gtf(annotation_dic):
    pass



def write2gff3(annotation_dic,is_evidence='auto',supress_defline_info=False):
    """returns a gff3 string following best gff3 standards (gene, mRNA, exon, and CDS lines) from an annotation dictionary"""
    outlines=[]
    for gene in annotation_dic:
        ghandle=annotation_dic[gene]
        outlines=outlines+annotation2gff3(ghandle,is_evidence=is_evidence,supress_defline_info=supress_defline_info).split('\n')
    return '\n'.join(outlines)

def write2appolo_gff3_scaffolds(annotation_dic,genomelocation='default',is_evidence='auto',supress_defline_info=False,split_scaffold_id=False):
    """returns a list of scaffold gff3 strings with fasta embedded"""
    annotation_dic.build_locus_dic()
    outlist=[]
    if genomelocation=='default':
        indi=list(annotation_dic)[0]
        try:
            genomefile=open(annotation_dic[indi].genomelocation)
        except:
            print "Cannot open genome file "+annotation_dic[indi].genomelocation
    else:
        genomefile=open(genomelocation)
    genome=genomefile.read().split('>')[1:]
    genomedic={}
    for i in genome:
        name=i.split('\n')[0]
        if split_scaffold_id:
            name=name.split()[0]
        seq=''.join(i.split('\n')[1:])
        genomedic[name]=seq
    for locus in annotation_dic.locus_dic:
        lhandle=annotation_dic.locus_dic[locus]
        scaflines=[]
        for geneid in lhandle:
            scaflines.append(annotation2gff3(annotation_dic[geneid],is_evidence=is_evidence,supress_defline_info=supress_defline_info))
        scaflines.append('#Fasta\n>'+locus+'\n'+genomedic[locus])
        outlist.append('\n'.join(scaflines))
    return outlist







#def longest_read_frame(nucl_sequence,output='nucl'):





