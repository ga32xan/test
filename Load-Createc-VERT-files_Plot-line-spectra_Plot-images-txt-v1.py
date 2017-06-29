################################
### Check scalebar size!
### Mathias Pörtner 
### Check for each line spectra tagged lines with "Adapt for each line spectra"
### written for python 3.6
################################
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pylab as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np
from scipy.optimize import curve_fit as cv
import glob
from matplotlib import gridspec
import re

def string_simplify(str):
    """simplifies a string (i.e. removes replaces space for "_", and makes it lowercase"""
    return str.replace(' ','_').lower()

def laden_spec(data):
	header = {}
	f = open(data, encoding='utf-8', errors='ignore')
	header_ended = False
	caption = re.compile(':*:')
	key = ''
	contents = ''
	while not header_ended:
		line = f.readline()
		if not line: break
		if line[0:4] == "    ":
			parts=line.split()
			posi=np.array([float(parts[-2]),float(parts[-1])],float)
			header_ended = True
		else:
			parts = line.split('=')
			if len(parts)!=2: continue
			key, contents = parts
			line = line.strip()
			key = string_simplify(key) # set new name
			header[key] = contents.strip() # [todo: add some parsing here
	f.close()
	
	dacstep=np.array([float(header['delta_x_/_delta_x_[dac]']),float(header['delta_y_/_delta_y_[dac]'])],float)
	pixelsize=np.array([float(header['num.x_/_num.x']),float(header['num.y_/_num.y'])],float)
	imagesize=np.array([float(header['length_x[a]']),float(header['length_y[a]'])],float)
	
	posi=posi/dacstep
	posi[0]=(pixelsize[0]/2.0+posi[0])*imagesize[0]/pixelsize[0]/10
	posi[1]=(pixelsize[1]-posi[1])*imagesize[1]/pixelsize[1]/10
    
	A=np.genfromtxt(data,delimiter='	',skip_header=212,skip_footer=0)
	U=A[:,3]
	dIdU=A[:,2]
	return(U,dIdU,posi)

def laden_image(data):
	header = {}
	f = open(data, encoding='utf-8', errors='ignore')
	header_ended = False
	caption = re.compile(':*:')
	key = ''
	contents = ''
	while not header_ended:
		line = f.readline()
		if not line: break
		if line[0] != "#":
			header_ended = True
		else:
			parts = line.split(':')
			if len(parts)!=2: continue
			key, contents = parts
			line = line.strip()
			key = string_simplify(key[2:])
			header[key] = contents[:-4].strip()
	f.close()
	
	ext=np.array([float(header['width']),float(header['height'])],float)
	
	X=np.loadtxt(data)*1e10
	
	return(X,ext)

def contrast(X,col):
	mini=100.0
	maxi=0.0
	for i in X:
		if max(i)>maxi:
			maxi=max(i)
		if min(i)<mini:
			mini=min(i)
	mean=np.mean([mini,maxi])

	global fig
	
	fig, ax = plt.subplots()
	plt.subplots_adjust(left=0.25, bottom=0.25)
	
	global unten
	global oben
	
	unten=mini
	oben=maxi

	global l
	
	l = plt.imshow(X, cmap=col, aspect='auto', vmin=mini, vmax=maxi)

	axcolor = 'lightgoldenrodyellow'
	axunten = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
	axoben = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
	
	global sunten
	global soben
	
	sunten = Slider(axunten, 'Unten', mini-(mean-mini), maxi, valinit=mini)
	soben = Slider(axoben, 'Oben', mini, maxi+(maxi-mean), valinit=maxi)
	sunten.on_changed(update)
	soben.on_changed(update)

	global buttonres
	
	resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
	buttonres = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

	global buttonsav
	
	saveax = plt.axes([0.65, 0.025, 0.1, 0.04])
	buttonsav = Button(saveax, 'Save', color=axcolor, hovercolor='0.975')
	
	buttonres.on_clicked(reset)
	buttonsav.on_clicked(save)
	buttonsav.on_clicked(change_color)

	plt.show()

	return(unten,oben)

def update(val):
	global sunten
	global soben
	global l
	global fig
	vunten = sunten.val
	voben = soben.val
	l.set_clim(vmin=vunten, vmax=voben)
	fig.canvas.draw_idle()

def reset(event):
	global sunten
	global soben
	sunten.reset()
	soben.reset()

def change_color(event):
	global buttonsav
	global fig
	buttonsav.color = 'red'
	buttonsav.hovercolor = buttonsav.color
	fig.canvas.draw()
    
def save(val):
	global unten
	global oben
	global sunten
	global soben
	unten = sunten.val
	oben = soben.val

def f(x,a,b,c):
	return(a*(x-b)**2+c)
	
#Load data
#Adapt for each line spectra
spec=glob.glob('F170509.100544.L*.VERT')
spec.sort()
#Adapt for each line spectra
ima,imagesize=laden_image('F170509.095922.txt')
#Adapt for each line spectra
didv,didvsize=laden_image('F170509.113636-didv.txt') #optional dI/dV channel

ext=[]
matrixx=[]
matrixy=[]
spec_posi=[]
for i in spec:
	x,y,posi=laden_spec(i)
	matrixy.append(y)
	matrixx.append(x)
	spec_posi.append(posi)
	if i==spec[0] or i==spec[-1]:
		ext.append(posi)
matrixx=np.array(matrixx,float)
matrixy=np.array(matrixy,float)

line_length=np.sqrt((ext[0][0]-ext[1][0])**2+(ext[0][1]-ext[1][1])**2)


#normalize specs
# change for range where max is found => -700 , 0
maxi=[]
globmaxy=0
for n,i in enumerate(matrixy):
	maxi.append(0)
	for m,j in enumerate(i):
	#Adapt for each line spectra
		if matrixx[n][m]>-700 and matrixx[n][m]<0 and matrixy[n][m]>maxi[n]:
			maxi[n]=matrixy[n][m]
	if max(i)>globmaxy:
		globmaxy=max(i)
#average specs, floating average for 5pt  so if m%5==4: should be m%(x)==(x-1):     and    [m-(x-1):m])/x
matrixyneu=[]
for n,i in enumerate(matrixy):
	matrixyneu.append([])
	for m,j in enumerate(i):
		if m%5==4:
			matrixyneu[-1].append(sum(matrixy[n][m-4:m])/5)
matrixyneu=np.array(matrixyneu,float)

# New or old contrast
#### contrast values are saved in *.csv
cons=[]
ans='-'
while ans!='o' or ans!='n':
	ans=input('Use old or new values for contrast? (o/n) ')
	if ans=='o':
	#Adapt for each line spectra
		cons=np.loadtxt('F170509.095922.csv',delimiter=',')
		break
	elif ans=='n':
		#find out clim topo
		topounten,topooben=contrast(ima,'viridis')
		cons.append(topounten)
		cons.append(topooben)
		#find out clim didv
		didvunten,didvoben=contrast(didv,'afmhot')
		cons.append(didvunten)
		cons.append(didvoben)
		#find out clim specs along line
		specunten,specoben=contrast(matrixyneu.T,'afmhot')
		cons.append(specunten)
		cons.append(specoben)
		#Adapt for each line spectra
		np.savetxt('F170509.095922.csv',cons,delimiter=',')
		break

#plot
fontna=16
fontnu=12
##ohne dI/dV => figsize=(10,5)
####figsize=(15,5) is three times wider than high, dpi
plt.figure(figsize=(15,5),dpi=100,tight_layout=True)
##ohne dI/dV => .GridSpec(1, 2)
gs = gridspec.GridSpec(1, 3)

##topo
plt.subplot(gs[0])
plt.title('Topography',fontsize=fontna)
plt.imshow(ima,cmap='viridis',extent=[0,imagesize[0],0,imagesize[1]],vmin=cons[0],vmax=cons[1])
##ohne achsenbeschriftung, mit achsenbeschritung folgende zwei zeilen auskommendtieren
plt.xticks([])
plt.yticks([])
#plots points of spectra in topography
for pos in spec_posi:
	plt.plot(pos[0],pos[1],'.',color='white',ms=1)
#sets scalebar (0,0) unten links referenz
	for x in np.linspace(1,6,1000):
	plt.plot(x,1,color='white')

##didv
plt.subplot(gs[1])
plt.title('dI/dV',fontsize=fontna)
plt.imshow(didv,cmap='afmhot',extent=[0,11.05,0,11.05*230/512],vmin=cons[2],vmax=cons[3])
##ohne achsenbeschriftung, mit achsenbeschritung folgende zwei zeilen auskommendtieren
plt.xticks([])
plt.yticks([])
#sets scalebar (0,0) unten links referenz
for x in np.linspace(1,6,1000):
	plt.plot(x,1,color='white')

##specs
plt.subplot(gs[2])
#configures line spectry 'map'
plt.imshow(matrixy.T,cmap='afmhot',extent=[0,line_length,min(matrixx[0]),max(matrixx[0])],aspect='auto',vmin=cons[4],vmax=cons[5])
#sets xlabel
plt.xlabel('Distance x [nm]',fontsize=fontna)
plt.xticks(fontsize=fontnu)
#sets y-label
plt.ylabel('Bias voltage [mV]',fontsize=fontna)
plt.yticks([-1000,0,1000],fontsize=fontnu)
#inserts color bar
plt.colorbar()
#saves fiogure in pdf file, png, jpg
plt.savefig("test.pdf")
#shows plot
plt.show()

