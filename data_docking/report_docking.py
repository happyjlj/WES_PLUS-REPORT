#!/share/public/software/Python-2.7.13/bin python2.7
#-*-coding=utf-8-*-
from optparse import OptionParser
#from ReportClass import *

import create_family_xml
import create_common_xml
import ConfigParser
import urllib2,urllib

import json  
import copy  
import glob 
import stat
import time
import sys
import os  
import re
global pdf_name
global absolute_json
global ISOTIMEFORMAT
ISOTIMEFORMAT='%Y-%m-%d %X'

absolute_json="/share/work4/jialj/REPORT_BACK_UP/json_online/"
pdf_name=""


##阴性附录疾病描述加粗
def get_bold_appendix(list):
	disease_bold=[]
	for item in list:
		if 'REMARK' in item:
			temp=item['REMARK'].split('。')
			for jtem in temp:
				disease=temp[0].replace('\n','').replace(' ','').split('的临床特征')[0]
				if disease not in disease_bold:
					disease_bold.append(disease)
	return disease_bold
###判断合并后的数据是否有NA的项，如果有则去掉,如果全都是NA，只保留第一个（默认是先证者）
def del_NA_in_table(list):
	tmp=[]
	if len(list)>=2:
		for item in list:
			if 'position' in item.keys() and  "NA" in item['position'] :
				tmp.append(item)
			else:
				pass
			if 'chr_position' in item.keys() and "NA" in item['chr_position']:
				tmp.append(item)
			else:
				pass	
	for i in tmp:
		list.remove(i)
	if len(list)==0:
		list.append(tmp[0])
	else:
		pass
	return list 



###获取家系其他成员的核心CNV信息
def get_other_family_members_cnv(sample_id,bus_code_combine,list):
	for id in  bus_code_combine:
		if id == sample_id:
			pass
		else:
			family_info={}
			url=r'http://10.100.11.55/berry/gdwes/getByBusCode/'+id
			try:
				html=urllib2.urlopen(url,timeout=10)
			except  urllib2.URLError,e:
				print "The clinical API did not return data, please check the url:'%s'!"%url
				sys.exit(1)
			sample_info=html.read()
			if sample_info:
				hjson = json.loads(sample_info)	#转为json数据
				result=hjson['success']
				if result==True:
					family_info = hjson
					if len(family_info['vus']):
						for item in family_info['vus'] :
							if 'NA' not in item['chr_position']:
								list.append(item)
	###判断合并后的数据是否有NA的项，如果有则去掉
	list=del_NA_in_table(list)
	return	list

###统一字段名称
def trans_ne_appendix_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			tmp['Sample']=item['sample']
			del(tmp['sample'])
			if 'gene' in item.keys():
				tmp['Gene']=item['gene']
				del(tmp['gene'])
			if 'position' in item.keys():
				tmp['Position']=item['position']
				del(tmp['position'])
			if 'NM' in item.keys():
				tmp['Transcript']=item['NM']
				del(tmp['NM'])
			if 'exon' in item.keys():
				tmp['Exon']=item['exon']
				del(tmp['exon'])
			if 'hgvs_c' in item.keys():
				tmp['hgvs.c']=item['hgvs_c']
				del(tmp['hgvs_c'])
			if 'hgvs_p' in item.keys():
				tmp['hgvs.p']=item['hgvs_p']
				del(tmp['hgvs_p'])
			if 'var_type' in item.keys():
				tmp['Type']=item['var_type']
				del(tmp['var_type'])
			##阴性附录的genotype特殊处理
			if 'genotype' in item.keys():
				tmp['Genotype']=item['genotype']
				del(tmp['genotype'])
			###其他家系的杂合性
			if 'father' in item.keys():
				if 'NA' == item['father']:
					tmp['Father']=''
				else:
					tmp['Father']=item['father']
				del(tmp['father'])
			else:
				tmp['Father']=''
			if 'mother' in item.keys():
				if 'NA' == item['mother']:
					tmp['Mother']=''
				else:
					tmp['Mother']=item['mother']
				del(tmp['mother'])
			else:
				tmp['Mother']=''
			if 'others' in item.keys():
				if 'NA' == item['others']:
					tmp['Others']=''
				else:
					tmp['Others']=item['others']
				del(tmp['others'])
			else:
				tmp['Others']=''
			################################
			if 'clinical_level' in item.keys():
				tmp['ACMGLevel']=item['clinical_level']
				del(tmp['clinical_level'])
			if 'disease' in item.keys():
				tmp['Disease']=item['disease']
				del(tmp['disease'])
			if 'inheritance' in item.keys():
				tmp['Inheritance']=item['inheritance']
				del(tmp['inheritance'])
			if 'maf' in item.keys():
				tmp['AlleleFrequency']=item['maf']
				del(tmp['maf'])
			if 'rs' in item.keys():
				tmp['dbSNP']=item['rs']
				del(tmp['rs'])
			if 'REMARK' in item.keys():
				tmp['REMARK']=item['REMARK']
				#del(tmp['rs'])
			result.append(tmp)
	
	return result				
#########将NA替换为NA
def trans_NA_to_ND(list):
	result=[]
	##core_data
	for item in list:
		for i in item.keys():
			if 'NA' == item[i]:
				item[i]='ND'
	result=list
	return result
####填充空的扩展信息
def fill_ND_to_Extend(sample_id):
	sites={}
	sites['Gene']="ND"
	sites['Sample']=sample_id
	sites['Type']="ND"
	sites['Position']="ND"
	sites['Transcript']="ND"
	sites['Exon']="ND"
	sites['hgvs.c']="ND"
	sites['hgvs.p']="ND"
	sites['Genotype']="ND"
	sites['dbSNP']="ND"
	sites['AlleleFrequency']="ND"
	sites['Disease']="ND"
	sites['Inheritance']="ND"
	sites['ACMGLevel']="ND"
	sites['Typical_age_of_onset']="ND"
	return sites
####删除散样扩展报告中与核心或提示重复的位点
def del_bulk_repeat_sites(list1,list2,list3,sample_id):
	sites_extend={}
	index=1
	extend_new=[]
	sites_core={}
	sites_note={}
	for item in list1:
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			#gene=item['Gene']
			pos=item['Position']
			tmp=pos
			aa=tmp.encode('unicode-escape').decode('string_escape')
			sites_core[aa]=1
	for jtem in list2:
		if 'ND' not in jtem['Gene'] and 'ND' not in jtem['Position']:
			#gene=jtem['Gene']
			pos=jtem['Position']
			tmp=pos
			bb=tmp.encode('unicode-escape').decode('string_escape')
			sites_note[bb]=1
	tmp_extend=[]
	for item in list3:
		index=0
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			#gene=item['Gene']
			pos=item['Position']
			tmp=pos
			if tmp in sites_core.keys() or tmp in sites_note.keys():
				pass
			else:
				tmp_extend.append(item)
		else:
			tmp_extend.append(item)
	##如果去重后的先证者信息为空了，需要用ND填充
	if len(tmp_extend)==0:
		sites=fill_ND_to_Extend(sample_id)
		extend_new.append(sites)
	else:
		for jtem in tmp_extend:
			extend_new.append(jtem)
	return extend_new


####删除家系扩展报告中与核心或提示重复的位点
def del_family_repeat_sites(list1,list2,list3,sample):
	sample_id=sample
	#relation=relationship
	sites_core={}
	sites_extend={}
	index=1
	extend_new=[]
	sites_note={}
	for item in list1:
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			gene=item['Gene']
			pos=item['Position']
			tmp=gene+pos
			aa=tmp.encode('unicode-escape').decode('string_escape')
			sites_core[aa]=1
	for jtem in list2:
		if 'ND' not in jtem['Gene'] and 'ND' not in jtem['Position']:
			gene=jtem['Gene']
			pos=jtem['Position']
			tmp=gene+pos
			bb=tmp.encode('unicode-escape').decode('string_escape')
			sites_note[bb]=1
	tmp_extend=[]
	for item in list3:
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			gene=item['Gene']
			pos=item['Position']
			tmp=pos
			if tmp in sites_core.keys() or tmp in sites_note.keys():
				pass
			else:
				tmp_extend.append(item)
		else:
			tmp_extend.append(item)
	return tmp_extend

	
#--------------获取科诺安cnv图片路径信息--------------------_#
def get_cnv_picture_result(cnv_picture_file=[],bus_code="",status="",config_path=os.path.split(os.path.realpath(__file__))[0]): 
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	cnv_from_it_prefix= config.get("sample_path", "CNV_FROM_IT_PATH")
	cnv_bit_path_prefix=config.get("sample_path", "CNV_BIT_PATH")
	cnv_to_it_prefix= config.get("sample_path", "CNV_TO_IT_PATH")
	#利用sftp获取CNV结果图
	if status.strip():
		#command="sftp -oPort=3033 "+ "wes@10.100.16.45:/zonghe/sharedisk/sharedisk/cnv/newImage/"+bus_code.strip()+"*/*.png  /share/work1/zhanglj/cnv_picture/"
		command="sftp -oPort=3033 "+ cnv_from_it_prefix + bus_code.strip()+"*/*.png "+cnv_bit_path_prefix
		print "CNV FROM IT:",command
		os.system(command)
		
	#path="/share/work1/zhanglj/cnv_picture/"+bus_code.strip()+"*"
	path=cnv_bit_path_prefix + bus_code.strip()+"*"
	file_list=glob.glob(path)
	if file_list:
		cnv_picture_file=file_list
		cnv_picture_file=file_list
		for item in cnv_picture_file:
			##传给IT
			#command="scp -P 3033 "+item +" gps@10.100.16.45:/zonghe/sharedisk/gps-shared/WES/cnv_picture/"
			command="scp -P 3033 "+item +" " +cnv_to_it_prefix
			print "SCP KNA:",command
			os.system(command)
	else:
		pass
	if not len(cnv_picture_file):
		cnv_picture_file=[]
	return cnv_picture_file 

	#返回拼接好之后的质控路径信息
def get_all_qc_path(bus_code_list,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	qc_prefix= config.get("sample_path", "QC_PATH")
	extend_prefix=config.get("sample_path","EXTEND_PATH")
	qc_path_list=[]
	result=""
	if len(bus_code_list):
		for bus_code in bus_code_list:
			str="*/analysis/"+bus_code+"/6.QC/QC/"+bus_code+".Core.gene.xls" ##修改文件，2018.12.14
			sample_path=os.path.join(qc_prefix,str)
			if len(glob.glob(sample_path))>0:
				qc_path_list.append(glob.glob(sample_path)[0])
			else:
				print "The QC file path: %s is not exists!"%(sample_path)
				sys.exit(1)
		if len(qc_path_list)==1:
			result=qc_path_list[0]
		elif len(qc_path_list) >1:
			result=",".join(qc_path_list)
		else:
			print "There is no qc path for any bus code!"
			sys.exit(1)
	else:
		print "There is no bus code information!"
		sys.exit(1)
	return result
###########获取扩展报告绝对路径
def get_absolute_extend_path(bus_code_list,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	extend_prefix=config.get("sample_path", "EXTEND_PATH")
	extend_path_list=[]
	result=""
	if len(bus_code_list):
		for bus_code in bus_code_list:
			str_extend="*/analysis/"+bus_code+"/5.Interpretation/ACMG/"+bus_code+"_Extended.xls"
			tmp_extend_path=os.path.join(extend_prefix,str_extend)
			##扩展绝对路径
			if len(glob.glob(tmp_extend_path))>0:
				if "_10G" in glob.glob(tmp_extend_path)[0]:
					extend_path_list.append(glob.glob(tmp_extend_path)[1])
				else:
					extend_path_list.append(glob.glob(tmp_extend_path)[0])
			else:
				print "The Extend file path: %s is not exists!"%(tmp_extend_path)
				sys.exit(1)
		if len(extend_path_list)==1:
			result=extend_path_list[0]
		elif len(extend_path_list) >1:
			result=",".join(extend_path_list)
		else:
			print "There is no extend path for any bus code!"
			sys.exit(1)
	else:
		print "There is no bus code information!"
		sys.exit(1)
	return result
###########获取qc的信息########
def get_all_qc_data(qc_path_string):
	qc_list=qc_path_string.split(",")
	
	result=[]
	for item in qc_list:
		sites={}
		sample_id=item.split("/")[-1].split(".Core.gene.xls")[0]
		sites["sample"] = sample_id
		sites['seq_project']=u'人类全外显子组'
		if os.path.exists(item):
			f=open(item)
			
			lines=f.readlines()
			for line in lines:
				if "PCT_TARGET_BASES_20X" in line:
					tmp_cover_20x=float(line.split()[-1].strip())*100
					tmp=format(tmp_cover_20x,'.2f')
					sites['cover_20x'] = str(tmp)+"%"
		else:
			print "The QC data %s is not exists"%(item)
			sys.exit(1)
		result.append(sites)
	return result
def deal_every_file_extend_site(site_line,sampleid):
	result=[]
	if len(site_line):
		sites={}
		line_data=site_line.split("\t")
		item=[data.strip() for data in line_data]
		sites['Gene']=item[0]
		sites['Sample']=sampleid
		sites['Type']=item[1]
		sites['Position']=item[2]
		sites['Transcript']=item[3]
		sites['Exon']=item[4]
		sites['hgvs.c']=item[5]
		sites['hgvs.p']=item[6]
		if item[7]=='het':
			sites['Genotype']=u'杂合'
		elif item[7]=='hom':
			sites['Genotype']=u'纯合'
		elif item[7]=='hem':
			sites['Genotype']=u'半合子'
		sites['dbSNP']=item[8]
		sites['AlleleFrequency']=item[9]
		sites['Disease']=item[10]
		sites['Inheritance']=item[11]
		sites['ACMGLevel']=item[12]
		sites['Typical_age_of_onset']=item[13]
	return sites
def get_all_extend_data(extend_path_string):
	extend_list=extend_path_string.split(",")
	sites=[]
	for item in extend_list:
		if os.path.exists(item):
			sampleid=item.split("_Extended.xls")[0].split("/")[-1]
			f=open(item)
			lines=f.readlines()
			sample_site={}
			if len(lines) >1:
				for line in lines[1:]:
					sites.append(deal_every_file_extend_site(line,sampleid))
			else:
				sample_site['Sample']=sampleid
				sample_site['Gene']='ND'
				sample_site['Type']='ND'
				sample_site['Position']='ND'
				sample_site['Transcript']='ND'
				sample_site['Exon']='ND'
				sample_site['hgvs.c']='ND'
				sample_site['hgvs.p']='ND'
				sample_site['Genotype']='ND'
				sample_site['dbSNP']='ND'
				sample_site['AlleleFrequency']='ND'
				sample_site['Disease']='ND'
				sample_site['Inheritance']='ND'
				sample_site['ACMGLevel']='ND'
				sample_site['Typical_age_of_onset']='ND'
				sites.append(sample_site)
		else:
			print "The Extend table %s is not exists"%(extend_file)
			sys.exit(1)
	return sites

###获取关键词
def get_key_words(sample_code,config_path=os.path.split(os.path.realpath(__file__))[0]):
	result=''
	key_words_list=[]
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	keywords_prefix= config.get("sample_path", "KEY_WORDS")
	str="wes_website/Phenolyzer/"+sample_code+"/cn_input.file" ##修改文件，2018.12.14
	keywords_path=os.path.join(keywords_prefix,str)
	if os.path.exists(keywords_path):
		f=open(keywords_path)
		lines=f.readlines()
		for line in lines:
			term=line.strip()
			if term not in key_words_list and term !='':
				key_words_list.append(term)
	result='、'.join(key_words_list)
	return result,key_words_list
##############获取GD接口的位点信息#############
def get_gd_info(sample_code,outfile,ana_type):
	url=r'http://10.100.11.55/berry/gdwes/getByBusCode/'+sample_code
	try:
		html=urllib2.urlopen(url,timeout=10)
	except  urllib2.URLError,e:
		print "The clinical API did not return data, please check the url:'%s'!"%url
		sys.exit(1)
	sample_info=html.read()
	exam_gd={}
	sample_list=[]
	exam_gd['qc_path']=[]
	exam_gd['qc_data']=[]
	exam_gd['extend_path']=[]
	exam_gd['extend_data']=[]
	exam_gd['articles']=[]
	exam_gd['conclusion']=[]
	exam_gd['verifyResult']=[]
	exam_gd['check_result']=''
	exam_gd['key_words']='' ##关键词，add20190604
	exam_gd['appendix_red']='' #add20190619
	exam_gd['appendix_bold']='' #add20190619
	input_info=""
	temp=""
	pattern="[-_a-zA-Z0-9]+"
	bus_code_combine=[]
	if sample_info:
		hjson = json.loads(sample_info)	#转为json数据
		result=hjson['success']
		if result==True:
			exam_gd = hjson
			####判断解读认领的检测类型是否与系统一致
			if len(exam_gd['clinic']):
				tmp_anatype=exam_gd['clinic'][0]['sub_test_items']
				if '单独先证者' in tmp_anatype:
					sys_anatype='bulk'
				if '夫妻' in tmp_anatype:
					sys_anatype='family'
				if '家系' in tmp_anatype:
					sys_anatype='family'
				if sys_anatype==ana_type:
					pass
				else:
					print '所认领的检测类型与系统不符，请确认后提交数据！'
					sys.exit(1)
	##------1 去除family_carry里的\n
	if len(exam_gd['result']):
		for item in exam_gd['result']:
			if '\n' in item['family_carry']:
				tmp=item['family_carry'].replace('\n','')
				item['family_carry']=tmp
			else:
				pass
	##------2 去掉note中的CNV信息
	tmp_dict=[]
	if len(exam_gd['note'])>=1:
		for item in exam_gd['note']:
			if item['variant_type']=='SNP':
				tmp_dict.append(item)
			else:
				pass
	exam_gd['note']=tmp_dict
	if len(exam_gd['note'])>=1:
		for item in exam_gd['note']:
			if item['variant_type']=='SNP' and '\n' in item['family_carry']:
				tmp=item['family_carry'].replace('\n','')
				item['family_carry']=tmp
			else:
				pass
	##-----3 增加家系成员的顺序，保证报告中的家系信息也是按照这个顺序
	exam_gd['relation_order']=[u'先证者',u'父亲',u'父親',u'母親',u'母亲',u'夫妻',u'丈夫',u'妻子',u'爷爷',u'奶奶',u'外公',u'姥爷',u'外婆',u'姥姥',u'女儿',u'儿子',u'哥哥',u'大哥',u'二哥',u'三哥',u'姐姐',u'大姐',u'二姐',u'三姐',u'弟弟',u'大弟',u'大弟弟',u'二弟',u'二弟弟',u'三弟',u'三弟弟',u'妹妹',u'大妹',u'大妹妹',u'二妹',u'二妹妹',u'三妹',u'三妹妹',u'小妹',u'姑母',u'姑姑',u'大姑',u'大姑姑',u'二姑',u'二姑姑',u'三姑',u'三姑姑',u'小姑',u'叔叔',u'大叔',u'二叔',u'三叔',u'小叔',u'伯伯',u'大伯',u'二伯',u'三伯',u'舅舅',u'大舅',u'二舅',u'三舅',u'小舅',u'姨',u'大姨',u'二姨',u'三姨',u'小姨',u'堂兄',u'堂弟',u'堂姐',u'堂妹',u'表兄',u'表弟',u'表姐',u'表妹',u'侄子',u'侄女',u'外甥',u'外甥女',u'孙子',u'孙女',u'外孙子',u'外孙女',u'胎儿1',u'胎儿2',u'弟弟或妹妹',u'姑父']
	order_relation_data=[]
	exam_gd['relationship']=""
	if len(hjson['clinic']):
		for item in hjson['clinic']:
			if item['bus_code']==sample_code:
				temp = item['family_test']
				sample=item
	else:
		print "The  Sample:%s is not exists in GD database, please check it!"%(sample_code)
		sys.exit(1)
	temp_relation=[]
	if temp:
		if "，" in temp:
			temp = temp.replace("，", ",")
		if "," in temp:
			temp_data = temp.split(",")
			####按照家系的顺序存取数据
			for key in exam_gd['relation_order']:
				for item in temp_data:
					if key == item.split(' ')[0]:
						order_relation_data.append(item)
			temp_data=order_relation_data
			if len(temp_data):
				temp_data = [item.strip() for item in temp_data]
				#####增加判断，客服现在录入家系时会有空格，如FU 8E1048FU0
				for item in temp_data:
					if ' ' in item:
						new_item=item.split(' ')[1]
						search=re.search(pattern,new_item)
						temp_sample=search.group()
						relation=item.split(' ')[0]
						temp_relation.append(relation)
					else:
						search=re.search(pattern,item)
						temp_sample=search.group()
					##去掉末尾带R0的编号，包括FUR0，MUR0,add20190521
					str_end=temp_sample.rfind('R0')
					if str_end!=-1:
						continue
					if ("V" in temp_sample) or ("MV" in temp_sample) or ("XV" in temp_sample) or ("XV1" in temp_sample) or ("XV2" in temp_sample) or ("XV3" in temp_sample) or ('-T' in temp_sample):
						continue
					bus_code_combine.append(temp_sample)
		else:
			search=re.search(pattern,temp)
			one_temp_sample=search.group()
			bus_code_combine=[one_temp_sample]
	##-----4 如果sample_code中有C0时，需要把family_test中的原编号忽略add20190521
	for item in bus_code_combine:
		if 'C0' in item and sample_code == item:
			tmp=item.replace('C0','')
			if tmp in bus_code_combine:
				bus_code_combine.remove(tmp)
		elif 'C0' in item:
			tmp=item.replace('C0','0')
			if tmp in bus_code_combine:
				bus_code_combine.remove(tmp)
	##-----5 判断合家系的位点的数据是否有NA的项，如果有则去掉
	exam_gd['result']=del_NA_in_table(exam_gd['result'])
	exam_gd['note']=del_NA_in_table(exam_gd['note'])
	##-----6 获取QC、Extend信息
	exam_gd['bus_code_combine'] = bus_code_combine
	qc_path_string=get_all_qc_path(bus_code_combine)
	extend_path_string=get_absolute_extend_path(bus_code_combine)
	exam_gd['qc_path'] = qc_path_string
	exam_gd['extend_path'] = extend_path_string
	exam_gd['qc_data'] = get_all_qc_data(qc_path_string)
	exam_gd['extend_data'] = get_all_extend_data(extend_path_string)
	##-----7 新增，关键词,附录关键词、疾病加粗，标红
	exam_gd['key_words'],exam_gd['appendix_red']=get_key_words(sample_code)
	exam_gd['appendix_bold']=get_bold_appendix(exam_gd['negative_appendix'])
	##-----8 获取subnumber
	for item in hjson['clinic']:
		if item['bus_code']==sample_code:
			exam_gd['subnumber']=item['subnumber']
	##----9 获取检测结论
	if hjson['result'][0]['result']=="negative":
		exam_gd['check_result']=u'阴性'
	elif hjson['result'][0]['result']=="positive":
		exam_gd['check_result']=u'阳性'
	elif hjson['result'][0]['result']=="unknown":
		exam_gd['check_result']=u'未知'
	else:
		exam_gd['check_result']=hjson['result'][0]['result']
	#-----10 WES+CNV时，提取科诺安图片地址
	sample_dict=[]
	find=False
	if len(exam_gd['clinic']):
		for item in exam_gd['clinic']:
			if item['bus_code']==sample_code.strip():
				sample_dict=item 
				find=True
	if find==False:
		print "The bus code:%s not found, please check the mongo url"%(sampleCode)
		sys.exit(1)
	exam_gd['imagepath']=sample_dict['imagepath'].strip()
	exam_gd['verify_cnv_data']=get_cnv_picture_result(bus_code=sample_code,status=exam_gd['imagepath'])
	exam_gd['verify_cnv_imgs']=trans_cnv_imgs(exam_gd['verify_cnv_data'])
	del exam_gd['verify_cnv_data']
	
	return exam_gd

def trans_vus_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			if 'chr_position' in item.keys():
				tmp['chr_position']=item['chr_position']
			if 'mut_size' in item.keys():
				tmp['mut_size']=item['mut_size']
			if 'disease' in item.keys():
				tmp['Disease']=item['disease']
				del(tmp['disease'])
			if 'sample' in item.keys():
				tmp['Sample']=item['sample']
				del(tmp['sample'])
			if 'mut_asses' in item.keys():
				tmp['mut_asses']=item['mut_asses']
			if 'result' in item.keys():
				tmp['result']=item['result']
			if 'bus_code' in item.keys():
				tmp['bus_code']=item['bus_code']
			if 'positive_reason' in item.keys():
				tmp['positive_reason']=item['positive_reason']
			if 'variant_type' in item.keys():
				tmp['variant_type']=item['variant_type']
			if 'type' in item.keys():
				tmp['type']=item['type']
			if 'gene_list' in item.keys():
				tmp['gene_list']=item['gene_list']
			result.append(tmp)
	
	return result
	
def trans_note_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			tmp['Sample']=item['sample']
			del(tmp['sample'])
			if 'gene' in item.keys():
				tmp['Gene']=item['gene']
				del(tmp['gene'])
			if 'position' in item.keys():
				tmp['Position']=item['position']
				del(tmp['position'])
			if 'NM' in item.keys():
				tmp['Transcript']=item['NM']
				del(tmp['NM'])
			if 'exon' in item.keys():
				tmp['Exon']=item['exon']
				del(tmp['exon'])
			if 'hgvs_c' in item.keys():
				tmp['hgvs.c']=item['hgvs_c']
				del(tmp['hgvs_c'])
			if 'hgvs_p' in item.keys():
				tmp['hgvs.p']=item['hgvs_p']
				del(tmp['hgvs_p'])
			if 'var_type' in item.keys():
				tmp['Type']=item['var_type']
				del(tmp['var_type'])
			if 'family_carry' in item.keys():
				tmp['Genotype']=item['family_carry']
				del(tmp['family_carry'])
				del(tmp['genotype'])
			if 'clinical_level' in item.keys():
				tmp['ACMGLevel']=item['clinical_level']
				del(tmp['clinical_level'])
			if 'disease' in item.keys():
				tmp['Disease']=item['disease']
				del(tmp['disease'])
			if 'inheritance' in item.keys():
				tmp['Inheritance']=item['inheritance']
				del(tmp['inheritance'])
			if 'maf' in item.keys():
				tmp['AlleleFrequency']=item['maf']
				del(tmp['maf'])
			if 'rs' in item.keys():
				tmp['dbSNP']=item['rs']
				del(tmp['rs'])
			result.append(tmp)
	
	return result
	
def trans_core_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			tmp['Sample']=item['sample']
			del(tmp['sample'])
			tmp['Gene']=item['gene']
			del(tmp['gene'])
			tmp['Position']=item['position']
			del(tmp['position'])
			tmp['Transcript']=item['NM']
			del(tmp['NM'])
			tmp['Exon']=item['exon']
			del(tmp['exon'])
			tmp['hgvs.c']=item['hgvs_c']
			del(tmp['hgvs_c'])
			tmp['hgvs.p']=item['hgvs_p']
			del(tmp['hgvs_p'])
			tmp['Type']=item['var_type']
			del(tmp['var_type'])
			tmp['Genotype']=item['family_carry']
			del(tmp['family_carry'])
			del(tmp['genotype'])
			tmp['ACMGLevel']=item['clinical_level']
			del(tmp['clinical_level'])
			tmp['Disease']=item['disease']
			del(tmp['disease'])
			tmp['Inheritance']=item['inheritance']
			del(tmp['inheritance'])
			tmp['AlleleFrequency']=item['maf']
			del(tmp['maf'])
			tmp['dbSNP']=item['rs']
			del(tmp['rs'])
			result.append(tmp)
	return result
def trans_verify_imgs(dict):
	index=1
	verify_imgs={}
	for item in dict:
		img_title=dict[item][1].split("/")[-1]
		verify_imgs[index]=[dict[item][0],img_title]
		index=index+1
	return verify_imgs
def trans_cnv_imgs(list):
	verify_cnv_img={}
	cnv_imgs={}
	imgs=[]
	for item in list:
		img_title=item.split("/")[-1]
		imgs.append(img_title)
	cnv_imgs=(',').join(imgs)
	return cnv_imgs

###提取真正的一代验证图片和wescnv图片
def get_verify_imgs_wescnv_imgs(dict,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	cnv_to_it_prefix= config.get("sample_path", "CNV_TO_IT_PATH")
	####
	wes_verify_imgs={}
	wescnv_imgs_names=''
	num=1
	verify_wescnv_imgs_list=[]
	##1)WES检测范围内的CNV
	tmp_wescnv_imgs_name=[]
	for item in dict.keys():
		tmp=[]
		title=dict[item][0]
		if 'WES检测范围' not in title:
			tmp=dict[item]
			wes_verify_imgs[num]=tmp
			num=num+1
		elif 'WES检测范围' in title:
			image=dict[item][1]
			###拷贝到IT服务器
			#command="scp -P 3033 " + image +" gps@10.100.16.45:/zonghe/sharedisk/gps-shared/WES/cnv_picture/"
			command="scp -P 3033 " + image +" "+cnv_to_it_prefix
			os.system(command)
			img_name=dict[item][1].split('/')[-1]
			tmp_wescnv_imgs_name.append(img_name)
			pass
		wescnv_imgs_names=(',').join(tmp_wescnv_imgs_name)
	return wes_verify_imgs,wescnv_imgs_names
	
def  main():
	usage = "Usage: %prog -i input_data_file -b busCode -a ana_type -c check_result -o output "
	parser=OptionParser(usage)
	parser.add_option("-i", "--input", dest="input", action="store", help="txt file contains the result information")
	parser.add_option("-b","--buscode",dest="buscode",action="store",help="The sample Business Code")
	parser.add_option("-a","--ana_type",dest="ana_type",help="The analysis type,option: bulk or family")
	parser.add_option("-c","--check_result",dest="check_result",help="The sample check result, option:positive or negative")
	parser.add_option("-o","--output",dest="output",help="The output directory name")
	parser.add_option("-y","--hospital",dest="hospital",default="jiahui",help="which hospital need the report, if the hospital is xinhua hospital,the parameter needed!")
	(options,args)=parser.parse_args()
	
	global pdf_name
	global CMD_path
	
	
	output_pdf_name=""
	all_json_file=""
	if not options.input:  
		print "The input data file contains the result information must be specified!"
		sys.exit(1)
	if not options.buscode:
		print "The business code of the sample must be specified!"
		sys.exit(2)
	if not options.ana_type or (options.ana_type!="bulk" and options.ana_type!="family"):
		print "The sample analysis type must be specified, options: bulk or family!"
		sys.exit(4)
	if not options.check_result or (options.check_result!="positive" and options.check_result!="negative"):
		print "The sample check result must be specified, options: negative or positive!"
		sys.exit(5)	
	if not options.output:
		output_pdf_name=options.buscode+".检测报告.pdf" 
		all_json_file=absolute_json+options.buscode+".all.json"###新增
	if options.output:
		if not os.path.exists(options.output):
			os.makedirs(options.output)
		output_pdf_name=options.buscode+".检测报告.pdf" 
		all_json_file=absolute_json+options.buscode+".all.json"###新增
	
	###1.命令调用输出
	CMD_path=os.path.abspath(sys.argv[0])
	input_path_prefix='/share/production/Genetics/WES/wes_website/Upload/*/'
	sample_input_path=input_path_prefix+options.buscode+".report.conf"
	abstract_path=glob.glob(sample_input_path)
	CMD='python2 '+CMD_path+' -i '+ ''.join(abstract_path) +' -b '+ options.buscode
	###2.分类模板#####################
	if options.ana_type=="bulk" and options.check_result=="positive":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、wescnv图片、文献等信息
		create_common_xml.buildNewsXmlFile(options.input)
		import read_xml
		supplemnt = read_xml.supplement_dictionary()
		#-------1.1 提取一代验证真正的图信息,以及WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		#----------1.1.1实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		#----2 获取gd 位点数据及qc，extend路径的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file,options.ana_type)
		##------2.1 获取输入模板的一代验证、wescnv、检测结论、方法、文献信息#############################################
		exam_gd['verify_result'] = trans_verify_imgs(wes_verify_imgs)
		exam_gd['verify_wescnv_imgs'] = wescnv_imgs_names
		exam_gd['articles'] = supplemnt['articles']
		exam_gd['conclusion_summary'] = supplemnt['conclusion_summary']
		exam_gd['check_conclusion'] = supplemnt['check_conclusion']
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['note_summary'] = supplemnt['note_summary']
		exam_gd['red'] = supplemnt['red']#add20190214
		exam_gd['overstriking'] =supplemnt['overstriking']#add20190214
		##---3 统一字段名称
		exam_gd['core_data']=trans_core_name(exam_gd['result'])
		del exam_gd['result']
		exam_gd['note_data']=trans_note_name(exam_gd['note'])
		del exam_gd['note']
		exam_gd['vus']=trans_vus_name(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		##---4 将NA改为ND
		exam_gd['core_data']=trans_NA_to_ND(exam_gd['core_data'])
		exam_gd['note_data']=trans_NA_to_ND(exam_gd['note_data'])
		exam_gd['vus']=trans_NA_to_ND(exam_gd['vus'])
		###add20190322
		exam_gd['negative_appendix']=trans_NA_to_ND(exam_gd['negative_appendix'])
		##---5 去重复,核心，提示语扩展报告中重复，去除扩展中的
		new_extend_data=del_bulk_repeat_sites(exam_gd['core_data'],exam_gd['note_data'],exam_gd['extend_data'],options.buscode)
		del exam_gd['extend_data']
		exam_gd['extend_data']=new_extend_data
		#print exam_gd['negative_appendix'][0].keys()
		##---6 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		##-----6.1 CMD调用
		if '阳性' in exam_gd['check_result']:
			CMD=CMD+' -a bulk -c '+'positive'
		print '对接程序：',CMD
		##---7 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		print "对接结束 : %s" % time.ctime()
	elif options.ana_type=="bulk" and options.check_result=="negative":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_common_xml.buildNewsXmlFile(options.input)
		import read_xml
		supplemnt = read_xml.supplement_dictionary()
		##---1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		##-------1.1.1 实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		##---2 获取gd 位点数据及qc，extend路径的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file,options.ana_type)
		##-----2.1 获取输入模板的检测结论、方法、文献信息#############################################
		exam_gd['verify_wescnv_imgs']=wescnv_imgs_names
		exam_gd['verify_result']=trans_verify_imgs(supplemnt['verify_result'])###转换一代验证
		exam_gd['articles']=supplemnt['articles']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']
		exam_gd['check_conclusion']=supplemnt['check_conclusion']
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['note_summary']=supplemnt['note_summary']
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214
		##---3 统一字段名称
		exam_gd['core_data']=trans_core_name(exam_gd['result'])
		del exam_gd['result']
		exam_gd['note_data']=trans_note_name(exam_gd['note'])
		del exam_gd['note']
		exam_gd['vus']=trans_vus_name(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		##---4 将NA改为ND
		exam_gd['core_data']=trans_NA_to_ND(exam_gd['core_data'])
		exam_gd['note_data']=trans_NA_to_ND(exam_gd['note_data'])
		exam_gd['vus']=trans_NA_to_ND(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_NA_to_ND(exam_gd['negative_appendix'])
		##---5 去重复,核心，提示语扩展报告中重复，去除扩展中的
		new_extend_data=del_bulk_repeat_sites(exam_gd['core_data'],exam_gd['note_data'],exam_gd['extend_data'],options.buscode)
		del exam_gd['extend_data']
		exam_gd['extend_data']=new_extend_data
		##---6 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		##----6.1 CMD调用
		if '阴性' in exam_gd['check_result']:
			CMD=CMD+' -a bulk -c '+'negative'
		if '未知' in exam_gd['check_result']:
			CMD=CMD+' -a bulk -c '+'unknown'
		print '对接程序：',CMD
		##---7 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		print "对接结束 : %s" % time.ctime()
	elif options.ana_type=="family" and options.check_result=="negative":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_family_xml.buildNewsXmlFile(options.input)
		import read_family_xml
		supplemnt=read_family_xml.supplement_dictionary()
		##---1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		##-------1.1.1 实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		##---2 获取gd 位点数据及qc，extend路径的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file,options.ana_type)
		##-----2.1 获取家系其他成员的CNV信息,合并
		exam_gd['vus']=get_other_family_members_cnv(options.buscode,exam_gd['bus_code_combine'],exam_gd['vus'])
		##-----2.2 获取输入模板的检测结论、方法、文献信息#############################################
		exam_gd['verify_wescnv_imgs']=wescnv_imgs_names
		exam_gd['verify_result']=trans_verify_imgs(supplemnt['verify_result'])###转换一代验证
		exam_gd['articles']=supplemnt['articles']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']
		exam_gd['check_conclusion']=supplemnt['check_conclusion']
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['note_summary']=supplemnt['note_summary']
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214
		##---3 统一字段名称
		exam_gd['core_data']=trans_core_name(exam_gd['result'])
		del exam_gd['result']
		exam_gd['note_data']=trans_note_name(exam_gd['note'])
		del exam_gd['note']
		exam_gd['vus']=trans_vus_name(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		##---4 将NA改为ND
		exam_gd['core_data']=trans_NA_to_ND(exam_gd['core_data'])
		exam_gd['note_data']=trans_NA_to_ND(exam_gd['note_data'])
		exam_gd['vus']=trans_NA_to_ND(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_NA_to_ND(exam_gd['negative_appendix'])
		##---6 去重复,核心，提示语扩展报告中重复，去除扩展中的
		new_extend_data=del_family_repeat_sites(exam_gd['core_data'],exam_gd['note_data'],exam_gd['extend_data'],options.buscode)
		del exam_gd['extend_data']
		exam_gd['extend_data']=new_extend_data
		##---7 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		##-----7.1 CMD调用#######
		if '阴性' in exam_gd['check_result']:
			CMD=CMD+' -a family -c '+'negative'
		if '未知' in exam_gd['check_result']:
			CMD=CMD+' -a family -c '+'unknown'
		print '对接程序：',CMD
		##---8 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		print "\n"
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		print "对接结束 : %s" % time.ctime()
		
	elif options.ana_type=="family" and options.check_result=="positive":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_family_xml.buildNewsXmlFile(options.input)
		import read_family_xml  
		supplemnt=read_family_xml.supplement_dictionary()
		##-----1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		##---------1.1.1 实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		##---2 获取gd 位点数据及qc，extend路径的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file,options.ana_type)
		##-----2.1 获取家系其他成员的CNV信息,合并
		exam_gd['vus']=get_other_family_members_cnv(options.buscode,exam_gd['bus_code_combine'],exam_gd['vus'])
		##-----2.2 获取输入模板的检测结论、方法、文献信息#############################################
		exam_gd['verify_wescnv_imgs']=wescnv_imgs_names
		exam_gd['verify_result']=trans_verify_imgs(supplemnt['verify_result'])###转换一代验证
		exam_gd['articles']=supplemnt['articles']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']
		exam_gd['check_conclusion']=supplemnt['check_conclusion']
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['note_summary']=supplemnt['note_summary']
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214
		exam_gd['appendix_bold']=get_bold_appendix(exam_gd['negative_appendix'])
		##---3 统一字段名称
		exam_gd['core_data']=trans_core_name(exam_gd['result'])
		del exam_gd['result']
		exam_gd['note_data']=trans_note_name(exam_gd['note'])
		del exam_gd['note']
		exam_gd['vus']=trans_vus_name(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		##---4 将NA改为ND
		exam_gd['core_data']=trans_NA_to_ND(exam_gd['core_data'])
		exam_gd['note_data']=trans_NA_to_ND(exam_gd['note_data'])
		exam_gd['vus']=trans_NA_to_ND(exam_gd['vus'])
		exam_gd['negative_appendix']=trans_NA_to_ND(exam_gd['negative_appendix'])
		##---6 去重复,核心，提示语扩展报告中重复，去除扩展中的
		new_extend_data=del_family_repeat_sites(exam_gd['core_data'],exam_gd['note_data'],exam_gd['extend_data'],options.buscode)
		del exam_gd['extend_data']
		exam_gd['extend_data']=new_extend_data
		##---7 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		##-----7.1 CMD#######
		if '阳性' in exam_gd['check_result']:
			CMD=CMD+' -a family -c '+'positive'
		print '对接程序', CMD
		##---8 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		print "对接结束 : %s" % time.ctime()
		
		
										
if	__name__=="__main__":
	main()   

