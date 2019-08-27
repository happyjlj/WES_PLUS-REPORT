#-*-coding:utf-8-*-
import xml.etree.ElementTree as ET
import sys
import os
import re
reload(sys)
sys.setdefaultencoding('utf8')

def indent(elem, level=0):
	i = "\n" + level*"  "
	if len(elem):
		if not elem.text or not elem.text.strip():
		  elem.text = i + "  "
		if not elem.tail or not elem.tail.strip():
		  elem.tail = i
		for elem in elem:
		  indent(elem, level+1)
		if not elem.tail or not elem.tail.strip():
		  elem.tail = i
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail = i
def buildNewsXmlFile(file):
	dict=read_family_info(file) 
	#建立根节点
	data = ET.Element("data")
	#g根节点下新建子节点
	#-------创建核心报告-------------#
	coreReport=ET.SubElement(data,"coreReport")
	createCoreReport(coreReport, dict['core'])
	#--------创建核心检测结论--------------#
	conclusion_summary=ET.SubElement(data,"conclusion_summary")
	conclusion_article_xml(conclusion_summary,dict['conclusion_summary'],"paragraph")
	#-----创建提示关注结论---------#
	note_summary=ET.SubElement(data,"note_summary")
	conclusion_article_xml(note_summary,dict['note_summary'],"paragraph")
	#-------CNV-seq检测结果---------#
	cnv_seq=ET.SubElement(data,"cnv_seq")
	conclusion_article_xml(cnv_seq,dict['cnv_seq'],"paragraph")
	
	#--------创建结论--------------#
	conclusion=ET.SubElement(data,"conclusion")
	conclusion_article_xml(conclusion,dict['conclusion'],"paragraph")
	#-------创建参考文献---------------#
	refArticle=ET.SubElement(data,"refArticle")
	#-------添加两篇固定的参考文献-------------#
	if not dict['refArticle']:
		dict['refArticle']=[]
	article_one="Kalia S S, Adelman K, Bale S J, et al. Recommendations for reporting of secondary findings in clinical exome and genome sequencing, 2016 update (ACMG SF v2.0): a policy statement of the American College of Medical Genetics and Genomics[J]. Genetics in Medicine Official Journal of the American College of Medical Genetics, 2017, 19(2):249."
	article_two="Richards S, Aziz N, Bale S, et al. Standards and Guidelines for the Interpretation of Sequence Variants: A Joint Consensus Recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology[J]. Genetics in Medicine Official Journal of the American College of Medical Genetics, 2015, 17(5):405."
	article_three=" Kearney H M, Thorland E C, Brown K K, et al. American College of Medical Genetics standards and guidelines for interpretation and reporting of postnatal constitutional copy number variants[J]. Genetics in Medicine Official Journal of the American College of Medical Genetics, 2011, 13(7):680-5."
	dict['refArticle'].insert(0,article_one)
	dict['refArticle'].insert(0,article_two)
	dict['refArticle'].append(article_three)
	conclusion_article_xml(refArticle,dict['refArticle'],"article")
	#-------------创建验证结果---------------#
	verifyResult=ET.SubElement(data,"verifyResult") 
	create_verify_xml(verifyResult,dict['verify'])
	indent(data)
	tree = ET.ElementTree(data)
	tree.write("common.xml", 'utf8')

#########创建验证结果的xml###############
def create_verify_xml(node,data):
	if len(data):
		for key in data.keys():
			sub = ET.SubElement(node, "result")
			verify_site=ET.SubElement(sub, "verify_site")
			verify_site.text=data[key][0]
			pic_name=ET.SubElement(sub, "pic_name")
			pic_name.text=data[key][1]
			
#####创建conclusion和article的xml##########
def conclusion_article_xml(node,data,label):
	if len(data):
		for item in data:
			sub=ET.SubElement(node,label)
			sub.text=item

###############创建扩展报告，label为添加的属性的key################
def createExtendReport(node,data):
	if len(data):
		index=0
		if len(data):
			for item in data:
				item_xml=ET.SubElement(node, "item")
				if len(item)==12:
					createTableXml(item_xml, item)


##############创建核心报告###################
def createCoreReport(node,data):
	if len(data):
		for item in data:
			subItem=ET.SubElement(node,"item")
			if len(item)==13:
				createTableXml(subItem, item,table="core")

########创建核心报告或扩展报告中的表的xml##################
def createTableXml(subItem,item,table="extend"):
	mut_gene = ET.SubElement(subItem, "mut_gene")
	mut_gene.text = "<font size='7'>"+item[0]+"</font>"
	mut_type = ET.SubElement(subItem, "mut_type")
	mut_type.text = "<font size='8'>"+item[1]+"</font>"
	mut_position = ET.SubElement(subItem, "mut_position")
	mut_position.text = "<font size='8'>"+item[2]+"</font>"
	transcript_code = ET.SubElement(subItem, "transcript_code")
	transcript_code.text = "<font size='8'>"+item[3]+"</font>"
	exon_code = ET.SubElement(subItem, "exon_code")
	exon_code.text = "<font size='8'>"+item[4]+"</font>"
	nucleotide_chagne = ET.SubElement(subItem, "nucleotide_chagne")
	nucleotide_chagne.text = "<font size='8'>"+item[5]+"</font>"
	acid_change = ET.SubElement(subItem, "acid_change")
	acid_change.text = "<font size='8'>"+item[6]+"</font>"
	hom_het = ET.SubElement(subItem, "hom_het")
	hom_het.text ="<font size='8'>"+ item[7]+"</font>"
	normal_frequency = ET.SubElement(subItem, "normal_frequency")
	normal_frequency.text = "<font size='8'>"+item[8]+"</font>"
	relate_disease = ET.SubElement(subItem, "relate_disease")	
	relate_disease.text = "<font size='8'>"+item[9]+"</font>"
	if item[9]=="ND":
		relate_disease.text = item[9]
	inherit_mode = ET.SubElement(subItem, "inherit_mode")
	inherit_mode.text = "<font size='8'>"+item[10]+"</font>"
	mut_evaluate = ET.SubElement(subItem, "mut_evaluate")
	mut_evaluate.text = "<font size='8'>"+item[11]+"</font>"  
	if table=="core":
		mut_origin=ET.SubElement(subItem, "mut_origin")
		mut_origin.text = "<font size='8'>"+item[12]+"</font>"	
def read_family_info(info_file):
	f=open(info_file,"r")
	content=f.read() 
	f.close()
	dict={}
	core=extract_file_content(content, pattern=r'<coreReport>(.*?)</coreReport>')
	conclusion=extract_file_content(content, pattern=r'<conclusion>(.*?)</conclusion>',default="conclusion")
	articles=extract_file_content(content, pattern=r'<article>(.*?)</article>',default="article")
	conclusion_summary=extract_file_content(content, pattern=r'<conclusion_summary>(.*?)</conclusion_summary>',default="conclusion")
	note_summary=extract_file_content(content, pattern=r'<note_summary>(.*?)</note_summary>',default="conclusion")
	cnv_seq=extract_file_content(content, pattern=r'<cnv_seq>(.*?)</cnv_seq>',default="conclusion")
	#extend=extract_extend_or_verify(content,default="extend")
	verify=extract_extend_or_verify(content,default="verify")
	dict['core']=core
	dict['conclusion_summary']=conclusion_summary
	dict['note_summary']=note_summary
	dict['conclusion']=conclusion
	dict['cnv_seq']=cnv_seq
	dict['refArticle']=articles
	#dict['extend']=extend
	dict['verify']=verify
	return dict
############提取读取文件内容中指定标签之内的内容################
def extract_file_content(content,pattern=r'<coreReport>(.*?)</coreReport>',default="core"):
	match=re.findall(pattern,content,re.M|re.S)
	if len(match):
		result=match[0]
		result=result.strip()
		if not result:
			result=""
			return result
		result_lines=result.split("\n")
		if default=="article":
			result=[data.strip() for data in result_lines]
		elif default=="conclusion":
			prefix=""
			res=[]
			num=0
			for para in result_lines:
				para=prefix+para
				temp=num  
				num+=len(para)
				#if(temp <2240 and num > 2240):
					#break_line=prefix+"<br/>"
					#res.append(break_line)
				res.append(para)
			result=res 
		elif default=="core": 
			res=[]
			for line in result_lines:
				line_data=line.split("\t")
				line_data=[data.strip() for data in line_data]
				res.append(line_data)
			result=res
	else:
		result=""
	return result

###########提取扩展报告和验证数据中的内容##########
def extract_extend_or_verify(content,default="extend",char_num=5):
	#------读取扩展报告中的内容-------#
	if default=="extend":
		pattern=r'<extendReport>(?P<item>.*?)</extendReport>'
		match = re.findall(pattern, content,re.S|re.M)
		if match:
			for item in match:
				all_sites=deal_extend_list(item)
		result=all_sites
		#print result
		#sys.exit(1)
	#-------读取验证数据内容------------#
	elif default=="verify":
		pattern = r'<result>(?P<item>.*?)</result>'
		match = re.findall(pattern, content, re.S | re.M)
		result={}
		if match:
			index=1
			for item in match:
				data_list= deal_verify_list(item, seperator="\n")
				if len(data_list) >1:
					result[index]=[data_list[0],data_list[1]]
					index+=1
		else:
			result=""
	else:
		result = ""

	return result

###########根据提取到的文本内容处理扩展报告数据#############
def deal_extend_list(string,seperator="\t"):
	string=string.strip()
	site_list=[]
	item_list=string.split("\n")
	item_list=[item.strip() for item in item_list]
	for item in item_list:
		data_list=item.split(seperator)
		data_list=[data.strip() for data in data_list]
		if len(data_list)!=12:
			print "Error:The extend table data should have 12 item, while you have "+str(len(data_list))
			sys.exit(1)
		site_list.append(data_list)
	return site_list

###########根据提取到的文本内容处理验证结果信息###########
def deal_verify_list(string,seperator="\n"):
	string=string.strip()
	data_list=[]
	if seperator in string:
		data_list=string.split("\n")  
		data_list=[data.strip() for data in data_list]
	return data_list

#buildNewsXmlFile("common.txt")
