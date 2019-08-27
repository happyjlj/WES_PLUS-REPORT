#-*-coding:utf-8-*-
import xml.etree.ElementTree as ET
import sys
import os
import re
import ConfigParser
reload(sys)
sys.setdefaultencoding('utf8')

try:
	tree = ET.parse("common.xml",parser=ET.XMLParser(encoding='utf-8'))		 #打开xml文档
	root = tree.getroot()				 #获得root节点
except Exception, e:
	print "Error:cannot parse file:common.xml."+str(e)
	sys.exit(1)

def get_core_report(root_node):
	core_node=root_node.find("coreReport")
	all_sites=[]
	for item in core_node.findall("item"):
		site=[]
		for child in item:
			site.append(child.text.strip())
		all_sites.append(site)
	return all_sites
def get_cnv_seq(root_node):
	con_node=root_node.find("cnv_seq")
	result=[]
	prefix=""
	for para in con_node:
		text=para.text
		if "®" in text:
			text=text.strip()
			#text=text.replace("®","<sup>®</sup>")
			text=prefix+text
		else:
			text=text.strip()
			text = prefix + text
		result.append(text)
	new_result="<para>".join(result)
	return new_result
def get_conclusion(root_node):
	con_node=root_node.find("conclusion")
	result=[]
	prefix=""
	for para in con_node:
		text=para.text
		if "®" in text:
			text=text.strip()
			#text=text.replace("®","<sup>®</sup>")
			text=text.replace("®"," ® ")
			text=prefix+text
		else:
			text=text.strip()
			text = prefix + text
		result.append(text)
	new_result="<para>".join(result)
	return new_result

def get_conclusion_summary(root_node):
	con_node=root_node.find("conclusion_summary")
	result=[]
	prefix=""
	for para in con_node:
		text=para.text
		if "®" in text:
			text=text.strip()
			#text=text.replace("®","<sup>®</sup>")
			text=prefix+text
		else:
			text=text.strip()
			text = prefix + text
		result.append(text)
	new_result="<para>".join(result)
	return new_result
def get_note_summary(root_node):
	con_node=root_node.find("note_summary")
	result=[]
	prefix=""
	for para in con_node:
		text=para.text
		if text==None:
			text=''
		else:
			if "®" in text:
				text=text.strip()
				#text=text.replace("®","<sup>®</sup>")
				text=prefix+text
			else:
				text=text.strip()
				text = prefix + text
		result.append(text)
	new_result="<para>".join(result)
	return new_result
		
def get_check_gene_list(root_node):
	con_node=root_node.find("check_gene_list")
	result=[]
	prefix=""
	for para in con_node:
		text=para.text
		text=text.strip()
		text = prefix + text
		result.append(text)
	new_result="<para>".join(result)
	return new_result
	
def get_articles(root_node):
	article_node=root_node.find("refArticle")
	result=""
	for article in article_node:
		if article.text:
			result+=article.text.strip()+"<para>"
	return result

def get_extend_report(root_node):
	extend_node=root_node.find("extendReport")
	result=[]
	for item in extend_node.findall("item"):
		site=[]
		for child in item:
			site.append(child.text.strip())
		result.append(site)
	return result

def get_verify_result(root,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	verify_to_it_prefix =config.get("sample_path", "VERIFY_TO_IT_PATH")
	root_node=root
	verify_node=root_node.find("verifyResult")
	result_dict={}
	index=1
	for result in verify_node:
		verify_site=result.find("verify_site").text.strip()
		verify_result=result.find("pic_name").text.strip()
		#拷贝给IT
		if 'WES检测范围' in verify_site:
			pass
		else:
			#command="scp -P 3033 " + verify_result +" gps@10.100.16.45:/zonghe/sharedisk/gps-shared/WES/wes_picture/"
			command="scp -P 3033 " + verify_result +" "+ verify_to_it_prefix
			
			print "SCP WES CNV:",command
			os.system(command)
		result_dict[index]=[verify_site,verify_result]
		index+=1
	return result_dict



def	supplement_dictionary():
	root_node=root
	supplement_result={}
	red=[]
	overstriking=[]
	##非贪婪匹配
	p1=re.compile(r'#(.*?)#')
	p2=re.compile(r'\<B\>(.*?)\<\/B\>')
	#core_data,已弃用
	supplement_result["core_data"]=get_core_report(root_node)
	#检测结果与解读
	supplement_result['conclusion_summary']=get_conclusion_summary(root_node)
	red1=p1.findall(supplement_result['conclusion_summary'])
	overstriking1=p2.findall(supplement_result['conclusion_summary'])
	#提示检测结果与解读
	supplement_result['note_summary']=get_note_summary(root_node)
	red2=p1.findall(supplement_result['note_summary'])
	overstriking2=p2.findall(supplement_result['note_summary'])
	#cnv结果
	supplement_result['cnv_seq']=get_cnv_seq(root_node)
	#文献
	supplement_result['articles']=get_articles(root_node)
	#检测结论
	supplement_result['check_conclusion']=get_conclusion(root_node)
	red3=p1.findall(supplement_result['check_conclusion'])
	#overstriking3=p2.findall(supplement_result['check_conclusion'])
	#标红合并、去重
	supplement_result['red']=list(set(red1+red2+red3))
	#加粗合并、去重
	supplement_result['overstriking']=list(set(overstriking1+overstriking2))
	##一代验证
	supplement_result['verify_result']=get_verify_result(root_node)
	return supplement_result

def main(file):
	try:
		tree = ET.parse(file, parser=ET.XMLParser(encoding='utf-8'))	# 打开xml文档
		root = tree.getroot()	# 获得root节点
	except Exception, e:
		print "Error:cannot parse file:%s, %s."%(file,str(e))
		sys.exit(1)
	return supplement_dictionary(root)

#get_verify_result()
