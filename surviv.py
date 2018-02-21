import re,os,sys,warnings,numpy,scipy,math,itertools;

from scipy import stats;
from numpy import *;
from multiprocessing import Pool;
from scipy.optimize import fmin_cobyla
from scipy.optimize import fmin_l_bfgs_b
from math import log;

numpy.random.seed(1231);
warnings.filterwarnings('ignore');

#obsolete variables: #ReadLength
read_length=1;
#JunctionLength
junction_length=1;
#splicing difference cutoff
cutoff=0.1;

#MultiProcessor
MultiProcessor=1;
if len(sys.argv)>=5:
	MultiProcessor=int(sys.argv[4]);

#survival time and events
time=[];event=[];Lambda=[];risk=[];
ifile=open(sys.argv[2]);
ifile.readline();ilines=ifile.readlines();
for i in ilines:
	element=re.findall('[^\t\n]+',i);
	if element[1]=='NA':
		time.append(-1);
		event.append(0);
	else:
		time.append(float(element[1]));
		event.append(float(element[2]));
for i in range(len(time)):
	d=sum((array(time)==time[i])*array(event));
	n=sum(array(time)>=time[i]);
	risk.append(float(d)/n);
#unique time point?
for i in range(len(time)):
	Lambda.append(sum((array(time)<=time[i])*array(risk)));

#update the baseline hazard 
def update_hazard(beta_cov,beta_psi,cov,psi,time,event,I,S,inc_length,skp_length):
	risk=[];Lambda=[];risk_unique=[0]*len(time);time_rec={};next_time_rec=[];
	cut=stats.mstats.mquantiles(time,array(range(11))/10.0);
	#print('cut');print(cut);print('beta_psi');print(beta_psi);
	for i in range(len(time)):
		adj=min(1,(sum(I[i])+sum(S[i]))*pow(10,-2));
		if sum(I[i])+sum(S[i])<=10:
			adj=adj*0.1;
		adj=1;
		temp=list(abs(time[i]-cut));index=temp.index(min(temp));
		if time[i]>cut[index]:
			this_time=cut[index];this_time_next=cut[index+1];
		else:
			this_time=cut[index-1];this_time_next=cut[index];
		next_time_rec.append(this_time_next);
		if this_time in time_rec:
			risk.append(time_rec[this_time]);continue;
		d=sum((array(time)>=this_time)*(array(time)<this_time_next)*array(event));
		risk_set=(array(time)>=this_time);
		#j is the index of the risk_set
		risk_set_beta=0;
		for j in range(len(risk_set)):
			if risk_set[j]==True:
				linear=exp(dot(array(beta_cov),array(cov[j]))+dot(array(beta_psi),array(veclogit(psi[j]))));
				weight=0;
				for e in range(len(I[j])):
					weight=weight+pow(beta_psi[e],2)/F(I[j][e],S[j][e],psi[j][e],inc_length,skp_length);
				weight=max(1,1/(1-weight))*adj;
				#weight=1-weight;
				#print('weight');print(weight);
				risk_set_beta=risk_set_beta+weight*linear;
		risk.append(float(d)/risk_set_beta);
		time_rec[this_time]=float(d)/risk_set_beta;
		risk_unique[i]=float(d)/risk_set_beta;
	#print('time_rec');print(time_rec);
	for i in range(len(time)):
		Lambda.append(sum((array(time)<=next_time_rec[i])*array(risk_unique)));
	return([risk,Lambda]);

#disabled in this version	
#clinical covariates
# if len(sys.argv)>=9:
	# ifile=open(sys.argv[8]);
	# ifile.readline();ilines=ifile.readlines();
	# cov=[];
	# for i in ilines:
		# element=re.findall('[^\t\n]+',i);
		# temp=[];
		# for j in element[1:]:
			# temp.append(float(j)+1);
		# cov.append(temp);
# else:
cov=[[]]*len(time);
#print('test_cov_input');print(cov);

#binomial MLE optimization functions
def logit(x):
	if x<0.001:
		x=0.001;
	if x>0.999:
		x=0.999;
	return(log(x/(1-x)));
	
def veclogit(x):
	res=[];
	for i in x:
		res.append(logit(i));
	return(res);
	
def F(I,S,psi,inc_length,skp_length):
	res=(2*psi-1-pow(psi,2))/pow(psi,2)/pow(1-psi,2);
	#res=0;
	res+=-1*I*skp_length*((2*inc_length+skp_length)*psi+skp_length*(1-psi))/pow(psi,2)/pow(inc_length*psi+skp_length*(1-psi),2);
	res+=-1*S*inc_length*((2*skp_length+inc_length)*(1-psi)+inc_length*psi)/pow(1-psi,2)/pow(inc_length*psi+skp_length*(1-psi),2);
	return(res);

#function to optimize the vector psi_k1, psi_k2 for all the exons of replicate k
def myfunc_l(x, *args):	
	I=args[0];S=args[1];cov=args[2];psi=args[3];risk=args[4];Lambda=args[5];event=args[6];
	effective_inclusion_length=args[7];effective_skipping_length=args[8];
	beta1_length=args[9];beta2_length=args[10];
	beta1=args[11][:beta1_length];beta2=args[11][beta1_length:(beta1_length+beta2_length)];
	psi_k=args[12];
	if sum(I[psi_k])+sum(S[psi_k])<=10:
		adj=1;
	else:
		adj=1;
	psi_replace=x;
	linear=dot(array(beta1),array(cov[psi_k]))+dot(array(beta2),array(psi_replace));
	#log-likelihood
	if risk[psi_k]==0:
		l1=event[psi_k]*(-1000)+event[psi_k]*linear-Lambda[psi_k]*exp(linear);
	else:
		l1=event[psi_k]*log(risk[psi_k])+event[psi_k]*linear-Lambda[psi_k]*exp(linear);
	l2=0;
	for e in range(len(beta2)):
		new_psi=effective_inclusion_length*psi_replace[e]/(effective_inclusion_length*psi_replace[e]+effective_skipping_length*(1-psi_replace[e]));
		l2+=I[psi_k][e]*log(new_psi)+S[psi_k][e]*log(1-new_psi);
	res1=l1+l2*adj;
	return(-1*res1);

#function to optimize the vector psi_k1, psi_k2 for all the exons of replicate k
def myfunc_l_der(x, *args):	
	I=args[0];S=args[1];cov=args[2];psi=args[3];risk=args[4];Lambda=args[5];event=args[6];
	effective_inclusion_length=args[7];effective_skipping_length=args[8];
	beta1_length=args[9];beta2_length=args[10];
	beta1=args[11][:beta1_length];beta2=args[11][beta1_length:(beta1_length+beta2_length)];
	psi_k=args[12];
	if sum(I[psi_k])+sum(S[psi_k])<=10:
		adj=1;
	else:
		adj=1;
	psi_replace=x;
	res1=[];
	linear=dot(array(beta1),array(cov[psi_k]))+dot(array(beta2),array(psi_replace));
	for e in range(len(beta2)):
		new_psi=effective_inclusion_length*psi_replace[e]/(effective_inclusion_length*psi_replace[e]+effective_skipping_length*(1-psi_replace[e]));
		new_psi_der=effective_inclusion_length*effective_skipping_length/pow(effective_inclusion_length*psi_replace[e]+effective_skipping_length*(1-psi_replace[e]),2);
		temp=event[psi_k]*beta2[e]-Lambda[psi_k]*beta2[e]*exp(linear);
		temp+=(I[psi_k][e]/new_psi*new_psi_der-S[psi_k][e]/(1-new_psi)*new_psi_der)*adj;
		res1.append(float(temp));
	return(-1*array(res1));
	
def myfunc_surv(x, *args):
	I=args[0];S=args[1];cov=args[2];psi=args[3];risk=args[4];Lambda=args[5];event=args[6];
	effective_inclusion_length=args[7];effective_skipping_length=args[8];
	beta1_length=args[9];beta2_length=args[10];output=args[11];
	beta1=x[:beta1_length];beta2=x[beta1_length:(beta1_length+beta2_length)];
	res1=0;res2=0;res1_l1=0;res1_l2=0;
	for i in range(len(event)):
		adj=min(1,(sum(I[i])+sum(S[i]))*pow(10,-2));
		if sum(I[i])+sum(S[i])<=10:
			adj=adj*0.1;
		adj=1;
		#print('test_myfunc_surv');print(beta1);print(cov[i]);print(beta2);print(psi[i]);
		linear=dot(array(beta1),array(cov[i]))+dot(array(beta2),array(veclogit(psi[i])));
		#log-likelihood
		if risk[i]==0:
			l1=event[i]*(-1000)+event[i]*linear-Lambda[i]*exp(linear);
		else:
			l1=event[i]*log(risk[i])+event[i]*linear-Lambda[i]*exp(linear);
		l2=0;
		for e in range(len(psi[i])):
			new_psi=effective_inclusion_length*psi[i][e]/(effective_inclusion_length*psi[i][e]+effective_skipping_length*(1-psi[i][e]));
			l2+=I[i][e]*log(new_psi)+S[i][e]*log(1-new_psi);
		#print('test_l2');print(l2);print(psi[i]);print(I[i][e]);print(S[i][e]);
		res1+=(l1+l2)*adj;
		#res1+=(l1+l2);
		res1_l1+=l1;res1_l2+=l2;
		#determinant
		F_prod=1;F_sum=0;
		for e in range(len(psi[i])):
			this_F=F(I[i][e],S[i][e],psi[i][e],effective_inclusion_length,effective_skipping_length);
			F_prod=F_prod*this_F;
			F_sum=F_sum+pow(beta2[e],2)/this_F;
		res2+=(-0.5*log(1-F_sum*Lambda[i]*exp(linear)))*adj;
		#res2+=(-0.5*log(1-F_sum*Lambda[i]*exp(linear)));
		#print('test_F_prod');print(F_prod);
	#if output==1:
	#	print('test_res');print(res1);print(res2);print(res1_l1);print(res1_l2);
	res1=float(res1);res2=float(res2);
	return(-1*(res2+res1));

def myfunc_surv_der(x, *args):
	I=args[0];S=args[1];cov=args[2];psi=args[3];risk=args[4];Lambda=args[5];event=args[6];
	effective_inclusion_length=args[7];effective_skipping_length=args[8];
	beta1_length=args[9];beta2_length=args[10];output=args[11];
	beta1=x[:beta1_length];beta2=x[beta1_length:(beta1_length+beta2_length)];
	#res1_deter: derivative of beta1 in the determinant; #res2_deter: derivative of beta2 in the determinant
	res1_deter=[];res2_deter=[];	
	#res1: derivative of beta1 in the log-likelihood; #res2
	res1=[];res2=[];
	for j in beta1:
		res1_deter.append(0);res1.append(0);
	for j in beta2:
		res2_deter.append(0);res2.append(0);
	#print('test_res');print(res1_deter);print(res2_deter);
	for i in range(len(event)):
		adj=min(1,(sum(I[i])+sum(S[i]))*pow(10,-2));
		if sum(I[i])+sum(S[i])<=10:
			adj=adj*0.1;
		adj=1;
		#res1_deter: derivative of beta1 in the determinant; #res2_deter: derivative of beta2 in the determinant
		this_res1_deter=[];this_res2_deter=[];
		linear=dot(array(beta1),array(cov[i]))+dot(array(beta2),array(veclogit(psi[i])));
		F_prod=1;F_array=[];
		for e in range(len(psi[i])):
			F_prod=F_prod*F(I[i][e],S[i][e],psi[i][e],effective_inclusion_length,effective_skipping_length);
			F_array.append(F(I[i][e],S[i][e],psi[i][e],effective_inclusion_length,effective_skipping_length));
		#determnant
		deter=(1-sum(pow(array(beta2),2)/array(F_array))*Lambda[i]*exp(linear));
		#beta2_der
		for e in range(len(beta2)):
			this_res2_deter.append((-1*Lambda[i])*exp(linear)*(2*beta2[e]/F_array[e]+sum(pow(array(beta2),2)/array(F_array))*logit(psi[i][e])));
			res2[e]+=(event[i]*logit(psi[i][e])-Lambda[i]*logit(psi[i][e])*exp(linear))*adj;
		#beta1_der
		for c in range(len(beta1)):
			this_res1_deter.append(sum(pow(array(beta2),2)/array(F_array))*(-1*Lambda[i])*exp(linear)*(cov[i][c]));
			res1[c]+=(event[i]*cov[i][c]-Lambda[i]*cov[i][c]*exp(linear))*adj;
		res1_deter=array(res1_deter)+(-0.5)*array(this_res1_deter)/deter*adj;
		res2_deter=array(res2_deter)+(-0.5)*array(this_res2_deter)/deter*adj;
	res=array(list(array(res1_deter)+array(res1))+list(array(res2_deter)+array(res2)))
	#print('test_res');print(res);print(res1_deter);print(res2_deter);
	return(-1*res);

def myfunc_surv_beta0(x, *args):
	I=args[0];S=args[1];cov=args[2];psi=args[3];risk=args[4];Lambda=args[5];event=args[6];
	effective_inclusion_length=args[7];effective_skipping_length=args[8];
	beta1_length=args[9];beta2_length=args[10];
	beta1=x[:beta1_length];beta2=array(x[beta1_length:(beta1_length+beta2_length)]);
	beta2[:]=0;
	res=myfunc_surv(list(beta1)+list(beta2),I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,beta1_length,beta2_length);
	return(res);

def myfunc_surv_beta0_der(x, *args):
	I=args[0];S=args[1];cov=args[2];psi=args[3];risk=args[4];Lambda=args[5];event=args[6];
	effective_inclusion_length=args[7];effective_skipping_length=args[8];
	beta1_length=args[9];beta2_length=args[10];
	beta1=x[:beta1_length];beta2=array(x[beta1_length:(beta1_length+beta2_length)]);
	beta2[:]=0;
	res=myfunc_surv_der(list(beta1)+list(beta2),I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,beta1_length,beta2_length);
	return(res[:beta1_length]);

def MLE_iteration(I,S,cov,time,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,beta_initial):
	psi=vec2psi(I,S,effective_inclusion_length,effective_skipping_length);
	beta_psi=[beta_initial]*len(I[0]);beta_cov=[0]*len(cov[0]);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;min_sum=pow(10,10);min_sum_beta_psi=0;
	#while(((abs(iter_cutoff)>0.001)&(abs(iter_cutoff)/(previous_sum+1)>pow(10,-6)))&(count<=iter_maxrun)):
	while (abs(iter_cutoff)>0.01)&(count<=iter_maxrun):
		if ((min_sum-previous_sum)<-2)&(count>10):
			break;
		#update baseline hazard
		res=update_hazard(beta_cov,beta_psi,cov,psi,time,event,I,S,effective_inclusion_length,effective_skipping_length);
		risk=res[0];Lambda=res[1];
		#print('test_risk');print(risk);print(sorted(Lambda));
		#update beta
		xopt=fmin_l_bfgs_b(myfunc_surv,list(beta_cov)+list(beta_psi),myfunc_surv_der,args=[I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,len(cov[0]),len(psi[0]),0],iprint=-1);
		#print('unconstrain_MLE_xopt');print(xopt);
		beta_cov=xopt[0][:len(cov[0])];beta_psi=xopt[0][len(cov[0]):];
		myfunc_surv(list(beta_cov)+list(beta_psi),I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,len(cov[0]),len(psi[0]),1);
		current_sum=xopt[1];
		if previous_sum!=0:
			iter_cutoff=current_sum-previous_sum;
		previous_sum=current_sum;count+=1;
		if min_sum>current_sum:
			min_sum=current_sum;
			min_sum_beta_psi=beta_psi;
		if sum(abs(array(beta_psi))>=20)>0:
			break;
		#update psi
		#if (iter_cutoff>0.1)&(count<=2):
		#	for i in range(len(event)):
		#		psi_init=[];psi_bound=[];
		#		for e in range(len(beta_psi)):
		#			psi_init.append(psi[i][e]);
		#			psi_bound.append([0.001,0.999]);
		#		xopt=fmin_l_bfgs_b(myfunc_l,psi_init,myfunc_l_der,args=[I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,len(cov[0]),len(psi[0]),list(beta_cov)+list(beta_psi),i],bounds=psi_bound,iprint=-1);
		#		psi[i]=list(xopt[0]);
		#	print('test_xopt');print(xopt);
		#print('test_psi');print(list(psi));
	#return([current_sum,[beta_cov,beta_psi,psi]]);
	return([min_sum,[beta_cov,min_sum_beta_psi,psi]]);

def MLE_iteration_constrain(I,S,cov,time,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,psi):
	psi=vec2psi(I,S,effective_inclusion_length,effective_skipping_length);
	beta_psi=[0]*len(I[0]);beta_cov=[0]*len(cov[0]);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;min_sum=pow(10,10);
	while (abs(iter_cutoff)>0.01)&(count<=iter_maxrun):
		#update baseline hazard
		res=update_hazard(beta_cov,beta_psi,cov,psi,time,event,I,S,effective_inclusion_length,effective_skipping_length);
		risk=res[0];Lambda=res[1];
		#print('test_risk');print(risk);print(sorted(Lambda));
		#update beta
		if len(cov[0])==0:
			xopt=myfunc_surv(beta_cov+beta_psi,I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,len(cov[0]),len(I[0]),1);
			current_sum=xopt;
		else:
			temp=[[-10,10]*(len(cov[0]))];
			xopt=fmin_l_bfgs_b(myfunc_surv_beta0,list(beta_cov),myfunc_surv_beta0_der,args=[I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,len(cov[0]),len(psi[0])],bounds=temp,iprint=-1);
			beta_cov=xopt[0][:len(cov[0])];
			current_sum=xopt[1];
		#print('constrain_MLE_xopt');print(xopt);
		if previous_sum!=0:
			iter_cutoff=current_sum-previous_sum;
		previous_sum=current_sum;count+=1;
		if min_sum>current_sum:
			min_sum=current_sum;
		#update psi
		#if (iter_cutoff>0.1)&(count<=2):
		#	for i in range(len(event)):
		#		psi_init=[];psi_bound=[];
		#		for e in range(len(beta_psi)):
		#			psi_init.append(psi[i][e]);
		#			psi_bound.append([0.05,0.95]);
		#		xopt=fmin_l_bfgs_b(myfunc_l,psi_init,myfunc_l_der,args=[I,S,cov,psi,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,len(cov[0]),len(psi[0]),list(beta_cov)+list(beta_psi),i],bounds=psi_bound,iprint=-1);
		#		psi[i]=list(xopt[0]);
				#print('constrain_MLE_xopt_psi');print(xopt);
		#print('test_psi');print(list(psi));
	return([min_sum,[beta_cov,beta_psi,psi]]);
	
#Random Sampling Function
def likelihood_test(i1,s1,cov,time,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,flag,id,beta_initial):
	#print('testing'+str(id));
	if flag==0:
		return([1,1]);
		#return([scipy.stats.uniform.rvs(0.1,0.9,1),1]);
	else:
		res=MLE_iteration(i1,s1,cov,time,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,beta_initial);
		res_constrain=MLE_iteration_constrain(i1,s1,cov,time,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,res[1][2]);
		#print('test');print(res);print(res_constrain);
		#return([1-scipy.stats.chi2.cdf(10*(abs(res_constrain[0]-res[0])),1)]);
		temp=scipy.stats.chi2.sf(2*(res_constrain[0]-res[0]),1);
		return([min(temp,1)]);

#MultiProcessorFunction
def MultiProcessorPool(n_original_diff):
	i1=n_original_diff[0];s1=n_original_diff[1];
	effective_inclusion_length=n_original_diff[2];effective_skipping_length=n_original_diff[3];
	time=n_original_diff[4];risk=n_original_diff[5];Lambda=n_original_diff[6];event=n_original_diff[7];
	cov=n_original_diff[8];
	flag=n_original_diff[9];id=n_original_diff[10];beta_initial=n_original_diff[11];
	P=likelihood_test(i1,s1,cov,time,risk,Lambda,event,effective_inclusion_length,effective_skipping_length,flag,id,beta_initial);
	return(P);

#Function for vector handling
def vec2float(vec):
	res=[];
	for i in vec:
		res.append(float(i));
	return(res);

#add 1 in both inclusion and skipping counts for robustness in small sample size
def vecAddOne(vec):
	res=[];
	for i in vec:
		res.append([i]);
	return(res);

def vecprod(vec):
	res=1;
	for i in vec:
		res=res*i;
	return(res);

def vecadd(vec1,vec2):
	res=[];
	for i in range(len(vec1)):
		res.append(vec1[i]+vec2[i]);
	return(res);
	
def vec2remove0psi(inc,skp):
	res1=[];res2=[];
	for i in range(len(inc)):
		if (inc[i]!=0) | (skp[i]!=0):
			res1.append(inc[i]);res2.append(skp[i]);
	return([res1,res2]);

def vec2psi_single(inc,skp,effective_inclusion_length,effective_skipping_length):
	psi=[];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	for i in range(len(inc)):
		temp=float(inc[i])/inclusion_length/(float(inc[i])/inclusion_length+float(skp[i])/skipping_length);
		if temp<=0.001:
			temp=0.001;
		if temp>0.999:
			temp=0.999;
		psi.append(temp);
	return(psi);	
	
def vec2psi(inc,skp,effective_inclusion_length,effective_skipping_length):
	psi=[];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	for i in range(len(inc)):
		psi.append([]);
		for e in range(len(inc[0])):
			temp=float(inc[i][e])/inclusion_length/(float(inc[i][e])/inclusion_length+float(skp[i][e])/skipping_length);
			if temp<=0.001:
				temp=0.001;
			if temp>0.999:
				temp=0.999;
			psi[i].append(temp);
	return(psi);


def vec210(vec):
	res=[];
	for i in vec:
		if i>0:
			res.append(1);
		else:
			res.append(-1);
	return(res);

def vec_remove_na_surv(inc,skp,time,risk,Lambda,event,cov):
	res_inc=[];res_skp=[];res_time=[];res_risk=[];res_Lambda=[];res_event=[];res_cov=[];
	for i in range(len(inc)):
		if (float(inc[i])+float(skp[i]))>0:
			if (time[i]>0):
				res_inc.append(inc[i]);
				res_skp.append(skp[i]);
				res_time.append(time[i]);
				res_risk.append(risk[i]);
				res_Lambda.append(Lambda[i]);
				res_event.append(event[i]);
				res_cov.append(cov[i]);
	return([res_inc,res_skp,res_time,res_risk,res_Lambda,res_event,res_cov]);
	
def myttest(vec1,vec2):
	if (len(vec1)==1) & (len(vec2)==1):
		res=stats.ttest_ind([vec1[0],vec1[0]],[vec2[0],vec2[0]]);
	else:
		res=stats.ttest_ind(vec1,vec2);
	return(res);

ifile=open(sys.argv[1]);
title=ifile.readline();
#analyze the title of the inputed data file to find the information of how much simulation are involved
#the min simulated round is 10, each time it increases by 10 times
element=re.findall('[^ \t\n]+',title);
ofile=open(sys.argv[3],'w');
ofile.write(title[:-1]+'\tPValue'+'\n');

list_n_original_diff=[];probability=[];psi_list_1=[];psi_list_2=[];
ilines=ifile.readlines();
for i in range(len(ilines)):
	element=re.findall('[^ \t\n]+',ilines[i]);
	if "NA" in ilines[i]:
		list_n_original_diff.append([0,0,0,0,0,0,0,0,0,0,element[0],0]);
		continue;
	inc=re.findall('[^,]+',element[1]);skp=re.findall('[^,]+',element[2]);
	effective_inclusion_length=int(element[3]);
	effective_skipping_length=int(element[4]);
	beta_initial=0;
	#if len(element)>=6:
	#	beta_initial=float(element[-1]);
	inc=vec2float(inc);skp=vec2float(skp);
	temp=vec_remove_na_surv(inc,skp,time,risk,Lambda,event,cov);
	inc_nona=temp[0];skp_nona=temp[1];
	psi_nona=vec2psi_single(inc_nona,skp_nona,effective_inclusion_length,effective_skipping_length)
	time_nona=temp[2];risk_nona=temp[3];Lambda_nona=temp[4];event_nona=temp[5];
	cov_nona=temp[6];
	#print('psi_nona');print(psi_nona);
	inc_95=stats.mstats.mquantiles(inc,0.95);
	skp_95=stats.mstats.mquantiles(skp,0.95);
#	if (((((len(inc_nona)<=20)|(sum(inc_nona)==0))|(sum(skp_nona)==0))|((len(set(psi_nona))/len(psi_nona))<=0.5))|(inc_95<=5))|(skp_95<=5):
	if (((((len(inc_nona)<=20)|(sum(inc_nona)==0))|(sum(skp_nona)==0))|((len(set(psi_nona)))<=10))|(inc_95<=3))|(skp_95<=3):
#		print('test');print(len(inc_nona));print(set(psi_nona));print(len(set(psi_nona))/len(psi_nona));print(inc_95);print(skp_95);
		list_n_original_diff.append([inc,skp,effective_inclusion_length,effective_skipping_length,time,risk,Lambda,event,cov,0,element[0],beta_initial]);
	else:
		inc=inc_nona;skp=skp_nona;
		inc=vecAddOne(inc);skp=vecAddOne(skp);
		list_n_original_diff.append([inc,skp,effective_inclusion_length,effective_skipping_length,time_nona,risk_nona,Lambda_nona,event_nona,cov_nona,1,element[0],beta_initial]);
	#if i>2:
	#	break;

if MultiProcessor>1:
	pool=Pool(processes=MultiProcessor);
	probability=pool.map(MultiProcessorPool,list_n_original_diff);
else:
	for i in range(len(list_n_original_diff)):
		probability.append(MultiProcessorPool(list_n_original_diff[i]));
index=0;
for i in range(len(ilines)):
    element=re.findall('[^ \t\n]+',ilines[i]);
    ofile.write(ilines[i][:-1]+'\t'+str(probability[i][0])+'\n');
ofile.close();
