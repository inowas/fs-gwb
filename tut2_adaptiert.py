
import numpy as np
import flopy

# Model domain and grid definition
Lx = 1000.
Ly = 1000.
ztop = 10.
zbot = -50.
nlay = 1
nrow = 20
ncol = 20
delr = Lx / ncol
delc = Ly / nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
hk = 1.
vka = 1.
sy = 0.1
ss = 1.e-4
laytyp = 1

# Variables for the BAS package
# Note that changes from the previous tutorial!
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
strt = 10. * np.ones((nlay, nrow, ncol), dtype=np.float32)

# Time step parameters
nper = 3
perlen = [1, 100, 100]
nstp = [1, 100, 100]
steady = [True, False, False]

# Flopy objects
modelname = 'tut2_adaptiert'
mf = flopy.modflow.Modflow(modelname, exe_name='mf2005')
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc,
                               top=ztop, botm=botm[1:],
                               nper=nper, perlen=perlen, nstp=nstp,
                               steady=steady)
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)
lpf = flopy.modflow.ModflowLpf(mf, hk=hk, vka=vka, sy=sy, ss=ss, laytyp=laytyp,
                               ipakcb=53)
pcg = flopy.modflow.ModflowPcg(mf)

# Make list for stress period 1,2,3
# Randbedingungen angepasst, linker Rand immer 10m , rechter Rand hier gerade 5m
# Eigentlich wollten wir das CHD package benutzen (weil wir das bisher für constant head immer benutzt habe,
# ...haben das aber nicht mit PloPy realisieren können
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
ibound[:, :, 0] = -1
ibound[:, :, -1] = -1
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
strt[:, :, 0] = 10.
strt[:, :, -1] = 5.
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

# Versuch CHD zu implementieren, hat aber nur funktioniert wenn wir ganze Matrix aufgeschrieben habe, ":" für "alle"
# ... funktioniert irgendwie nicht - für unterschiedliche Diskretisierung ist die Änderung dann nicht möglich/nervig
# stress_period_data ={0: [
  #  [0, ':', 0, 10., 10.], #west
   # [0, ':', -1, 9., 5.], #ost

    #],
# }
# chd = flopy.modflow.ModflowChd(mf, stress_period_data=stress_period_data)

# recharge konstant über gesamte Zeit
rch = flopy.modflow.ModflowRch(mf, nrchop=3, rech=1.2e-4) # Menge gerade noch willkürlich gewählt

# recharge pro Periode definiert?!?!?!

# Create the well package
# Remember to use zero-based layer, row, column indices!
# zweiten Brunnen implementiert (im oberen Viertel)
pumping_rate1 = -500.
pumping_rate2 = -300.
stress_period_data = {0: [[0, nrow/2 - 1, ncol/2 - 1, 0],
                          [0, nrow/4 - 1, ncol/2 - 1, 0]],
                      1: [[0, nrow / 2 - 1, ncol / 2 - 1, pumping_rate1],
                          [0, nrow / 4 - 1, ncol / 2 - 1, pumping_rate2]],
                      2: [[0, nrow / 2 - 1, ncol / 2 - 1, 0],
                          [0, nrow / 4 - 1, ncol / 2 - 1, 0]]}
wel = flopy.modflow.ModflowWel(mf, stress_period_data=stress_period_data)

# Output control
# Die folgenden Abschnitte sind eigentlich komplett aus dem Tutorial_2 entnommen, da wir die Befehle z.T. nicht verstehen
stress_period_data = {(0, 0): ['save head',
                               'save drawdown',
                               'save budget',
                               'print head',
                               'print budget']}
save_head_every = 1
oc = flopy.modflow.ModflowOc(mf, stress_period_data=stress_period_data,
                             compact=True)
# Wir würden gerne den zeitlichen Verlauf der heads an bestimmten Beobachtungspunkten (z.B. an beiden Brunnen) ausgeben - Wie?!
# Benutzen würden wir das HOB package, das einfügen funktioniert aber nicht bzw. bekommen wir dann immer Programmfehler
# Könnten Sie uns da weiterhelfen?!

# Write the model input files
mf.write_input()

# Run the model
success, mfoutput = mf.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')


# Imports
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

# Create the headfile and budget file objects
headobj = bf.HeadFile(modelname+'.hds')
times = headobj.get_times()
cbb = bf.CellBudgetFile(modelname+'.cbc')

# Setup contour parameters
levels = np.linspace(0, 10, 11)
extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
print('Levels: ', levels)
print('Extent: ', extent)

# Well point
# Diese sitzen nie genau über den Brunnen (zumindest in der Ausgabe die wir von FloPy bekommen)
# Was machen wir falsch?
wpt1 = ((float(ncol/2)-0.5)*delr, (float(nrow/2-1)+0.5)*delc)
wpt1 = (450., 550.)
wpt2 = ((float(ncol/4)-0.5)*delr, (float(nrow/2-1)+0.5)*delc)
wpt2 = (450., 750.)

# Make the plots
mytimes = [1.0, 101.0, 201.0]
for iplot, time in enumerate(mytimes):
    print('*****Processing time: ', time)
    head = headobj.get_data(totim=time)
    #Print statistics
    print('Head statistics')
    print('  min: ', head.min())
    print('  max: ', head.max())
    print('  std: ', head.std())

    # Extract flow right face and flow front face
    frf = cbb.get_data(text='FLOW RIGHT FACE', totim=time)[0]
    fff = cbb.get_data(text='FLOW FRONT FACE', totim=time)[0]

    #Create the plot
    plt.subplot(1, len(mytimes), iplot + 1, aspect='equal')
    plt.subplot(1, 1, 1, aspect='equal')
    plt.title('stress period ' + str(iplot + 1))


    modelmap = flopy.plot.ModelMap(model=mf, layer=0)
    qm = modelmap.plot_ibound()
    lc = modelmap.plot_grid()
    cs = modelmap.contour_array(head, levels=levels)
    plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
    quiver = modelmap.plot_discharge(frf, fff, head=head)


    mfc = 'None'
    if (iplot+1) == len(mytimes):
        mfc='black'
    plt.plot(wpt1[0], wpt1[1], lw=0, marker='o', markersize=8,
             markeredgewidth=0.5,
             markeredgecolor='black', markerfacecolor=mfc, zorder=9)
    plt.text(wpt1[0]+25, wpt1[1]-25, 'well', size=12, zorder=12)
    plt.plot(wpt2[0], wpt2[1], lw=0, marker='o', markersize=8,
             markeredgewidth=0.5,
             markeredgecolor='black', markerfacecolor=mfc, zorder=9)
    plt.text(wpt2[0] + 25, wpt2[1] - 25, 'wel2', size=12, zorder=12)
    plt.show()

plt.show()

# Plot the head versus time - wel1
idx = (0, nrow/2 - 1, ncol/2 - 1)
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Head at cell ({0},{1},{2})'.format(idx[0] + 1, idx[1] + 1, idx[2] + 1)
plt.title(ttl)
plt.xlabel('time')
plt.ylabel('head')
plt.plot(ts[:, 0], ts[:, 1])
plt.show()
# Plot the head versus time - wel2
idx = (0, nrow/4 - 1, ncol/2 - 1)
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Head at cell ({0},{1},{2})'.format(idx[0] + 1, idx[1] + 1, idx[2] + 1)
plt.title(ttl)
plt.xlabel('time')
plt.ylabel('head')
plt.plot(ts[:, 0], ts[:, 1])
plt.show()
# wenn wel1 oder wel2 direkt auf einer Kante zwischen 2 Zellen liegt bekommen wir eine Ausgabe von 0 über den ganzen Zeitraum
# ...Fehler von uns oder kann das Programm die Heads nicht berechnen, wenn der wel direkt auf einer Zellkante liegt?