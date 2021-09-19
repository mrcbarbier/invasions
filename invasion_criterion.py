# -*- coding: utf-8 -*-
from imports import *
from analysis import *
from cavity import *
from graphs import *
from models import *
from datatools import *
import pandas as pd
import matplotlib.figure as mpfig


Model=DynamicModel

def foodshow(meps):
    from reduction import measure_groups,plot_trophic_graph
    plot(meps.results['n'].matrix[:,:-1],hold=1,log='xy')
    plot(meps.results['n'].matrix[:,-1],hold=1,log='xy',lw=2,color='k')
    model,measure=meps,{}
    measure.update(model.export_params())
    measure_trophic(model,measure,
        groupby='trophic')
    plot_trophic_graph(model,measure,hold=1,nsize=0.04,minsize=0.0000001,niche='trophic')
    plt.show()


dft_dyn={
    'nlv':{'type':'simplelv',
        'variables':[('n','com'),],
        },
    }

dft_prm={
        'n':{
    'type':'variable',
    'axes': ('com',),
    'death': 10.**-8,
    'shape':(10,),
    'mean':1.,
    'std':1.,
    'sign':1,
    },

        'growth':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'growth',
    'mean':1,
    'std':0.2,
    },
        'selfint':{
    'type':'matrix',
    'variables':['n'],
    'dynamics':'nlv',
    'role':'diagonal',
    'mean':1,
    'std':0,
    },

        'community':{
    'type':'matrix',
    'variables':[('n','com'),('n','com') ],
    'role':'interactions',
    'mean':-.5,
    'stdrel':0.25,
    'dynamics':'nlv',
    'symmetry':1.,
    'diagonal':0,
    },

}
def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a



def invasibility(path='',res=60,rerun=1,extend=0,ncomm=1,ninv=1000,tmax=100000,tsample=1000,fname='invasibility',
                 mode='triangle',S=60,
                 gen_invaders=1, #Generate new invaders instead of using extinct species
                 **kwargs):
    path=Path(path)
    path.mkdir()

    prm=deepcopy(dft_prm)
    prm.update(kwargs.get('prm',{}))
    dyn=deepcopy(dft_dyn)
    dyn.update(kwargs.get('dyn',{}))

    if mode=='basic':
        gammas=np.linspace(-1,1,res)
        mus=np.ones(res)*.5
        sigmas=np.ones(res) *.1
        xs=gammas
    elif mode =='triangle':
        mus= np.tile(np.linspace(0,1,res+2)[1:-1],(res,1))
        delta=np.array([ np.linspace(0,min(mu,1-mu),res+1)[1:]  for mu in mus[0]]).T
        sigmas=delta/np.sqrt(3)
        mus=mus.T.ravel()
        sigmas=sigmas.T.ravel()
        gammas=np.ones(mus.shape)
        xs=zip(mus,sigmas)
    elif 'sigmadiff' in mode or 'siggam' in mode:

        prm['n']['shape']=(S,)
        res=20
        sigmas=np.linspace(0.05,.8,res)
        gammas=np.array([-.99,0,.99])
        mus = np.array([1,30])

        def cross(args):
            for i in range(len(args)):
                lst=[np.ones(a.shape ) for a in args[:i]]
                lst+=[args[i]]
                lst+=[np.ones(a.shape) for a in args[i+1:]]
                res=np.array([1])
                for l in lst:
                    res=np.multiply.outer(res,l)
                yield res.squeeze().ravel()

        sigmas,gammas,mus=cross([sigmas,gammas,mus])
    elif 'evol' in mode:
        """DEMONSTRATION MODE FOR ECOSYSTEM EVOLUTION WORKSHOP - PLOT OF W,V FIGURE"""
        prm['n']['shape']=(S,)
        res=4
        if 0:
            mus= np.tile(np.linspace(0,.2,res+2)[1:-1],(res,1))
            delta=np.array([ np.linspace(0,min(mu,1-mu),res+1)[1:]  for mu in mus[0]]).T
            sigmas=delta/np.sqrt(3)
            mus=mus.T.ravel()
            sigmas=sigmas.T.ravel()
        elif 'rndinv' in mode:
            mus = np.ones(1)*5.
            sigmas=np.ones(1)*.8
        else:
            mus = np.ones(1)*1.2
            sigmas=np.ones(1)*.88
        gammas=np.ones(mus.shape)*.99
        if 'foodweb' in mode:
            gammas=np.ones(mus.shape)*.0

        xs=zip(mus,sigmas)
        fname=mode #'gdrevoltab'
        if 'resource' in mode:
            for i in ('community','growth','selfint'):
                prm[i]['structure']='resource'
                prm[i].setdefault('requires',[])
                prm[i]['requires']+=['consumption']
            prm['growth']['requires']+=['r','mortality']
            prm['mortality']={
                'type':'matrix',
                'variables':['n'],
                'role':'mortality',
                'mean':0.1, 'std':0,
                }
            prm['consumption']={
                'type':'matrix',
                'variables':[('n','com'),('r','com') ],
                'mean':1., 'std':1./3.,
                'sign':1,
                'distribution':'exponential',
                'role':'consumption',
                'structure':'resource',
                'requires':['r'],
                }
            prm['r']={
                'type':'constant',
                'axes': ('com',),
                'role':'resource',
                'shape':(40,),
                'mean':1.,
                'std':0.,
                'sign':1,
                #'scaling':{'std':-.5},
                }


    table=[]
    PRMIDX=0

    if extend:
        try:
            table+=[row.to_dict() for i, row in pd.read_csv(path + fname + '.csv').iterrows()]
            print 'EXTENDED',len(table)
        except:
            pass
    elif not rerun:
        try:
            table = pd.read_csv(path + fname + '.csv')
            mus,sigmas,gammas=[],[],[]
        except:
            pass
    for mu,sigma,gamma in zip(mus,sigmas,gammas):
        PRMIDX+=1
        print 'Mu,sigma,gamma',mu,sigma,gamma
        for com in range(ncomm):
            print "Community",com
            S=prm['n']['shape'][0]
            if not 'foodweb' in mode:
                prm['community']['mean']=-mu/S
                if 'stdrel' in prm['community']:
                    del prm['community']['stdrel']
                prm['community']['std']=sigma/np.sqrt(S)
                # print '###',prm['community']['mean'],prm['community']['std']
                prm['community']['symmetry']=gamma
                prm['community']['distribution']='uniform'
            dpath=Path(path+'comm_{}'.format(com) )

            m=Model(parameters=prm,dynamics=dyn,**kwargs)
            #m.set_init(n=x0)
            tsample=10.
            m.evol(tmax=tmax,tsample=tsample,converge=1,reseed=False)
            # print m.data['niches'][m.results['n'].matrix[-1]>10**-5]

            # code_debugger()
            m.save(dpath)
            death = m.parameters['n']['death']

            measure={}
            measure_gen(m,measure,typ=['usual','effective'])
            measure.update(m.export_params())
            #measure_bunin(m,measure)
            # plot(m.results['n'].matrix,newfig=1)
            Nf=m.results['n'][-1]
            alive=np.where(Nf>death)[0]
            Nalive=Nf[alive]
            S=len(Nf)
            Salive=len(alive)
            if 'nonlin' in mode:
                refx = m.get_labeled_result(idx=-1)
                A=-np.dot(np.diag(1/Nf),m.get_jac(0,refx['n'].matrix ))
                np.fill_diagonal(A,0)
            else:
                A=-m.data['community'].matrix
            # code_debugger()
            Al=A[np.ix_(alive,alive)]
            ralive=m.data['growth'].matrix[alive]

            I=np.eye(len(alive))

            from scipy.linalg import inv as lainv,norm as lanorm, eig as laeig
            dNdK=lainv(I+ Al )
            worst=lanorm(dNdK,ord=2)
            vtrace, vsum = np.trace(dNdK), np.sum(dNdK)
            
            eigval,eigleft,eigright=laeig(np.dot(dNdK.T,dNdK),left=True,right=True)
            largesteig=np.argmax(np.abs(eigval))
            worstAi0=np.real(eigright[:,largesteig])
            worstA0i=np.dot(dNdK,worstAi0)/np.sqrt(np.real(eigval[largesteig])) #eigleft[:,largesteig]
            inv=0
            if 'obs' in mode:
                # observed=sorted(np.random.choice(list(range(Al.shape[0])),size=10))
                observed=sorted(np.argsort(np.argsort(Nalive))[-10:] )
                notobs=[z for z in range(Al.shape[0]) if not z in observed]
                Vnotobs=lainv(I[np.ix_(notobs,notobs)] + Al[np.ix_(notobs,notobs)] )
            else:
                observed=list(range(Al.shape[0]))
                notobs=[]
                Vnotobs=1.
            if  not gen_invaders:
                #use extinct species as invaders!
                ninv=S-Salive
            while inv<ninv:
                #print "Invader",inv
                if gen_invaders:
                    mgrowth=m.data['growth'].matrix
                    ginv=None
                    if 'foodweb' in mode and 'nontrophic' in mode:
                        r0=np.random.choice(mgrowth[mgrowth>0])*2.
                        ginv=0.
                    elif 'foodweb' in mode:
                        niches=m.data['niches'].matrix
                        niche0=np.random.choice(niches)
                        r0=np.random.choice(mgrowth[niches==niche0])
                        ginv=.7
                    else:
                        r0=np.random.normal(1,kwargs.get('invader_stdK', np.std(mgrowth) ))
                    if 'rndinv' in mode:
                        # ginv=np.random.uniform(-1,1)
                        if ginv is None:
                            ginv=1.
                        mean= (r0/ np.sum(Nalive)    ) * np.random.uniform(.7,1.3) #10**(np.random.uniform(-.5,.5))
                        # stdrel=np.random.uniform(0.01,.5)
                        # std=stdrel*mean
                        # std =np.sqrt( (mean**2 *vsum)/(vtrace)  ) * np.random.uniform(0,2.)
                        std=np.sqrt(1./vtrace *max(.5,1- mean**2 * vsum)  )* np.random.uniform(0,2.)
                        if 'nonlin' in mode:
                            std=std*2.  # I DONT KNOW WHY NORMAL STD MAKES SMALL VVARIATION
                            mean=mean

                        print 'Mean, mu/S, std,sigma/S**0.5',mean,mu/S, std,sigma/np.sqrt(S)
                    elif 'sigmadiff' in mode:
                        mean=-kwargs.get('invader_mean', m.get_param('community_mean') )
                        ginv=np.random.choice(gammas)
                        std=np.random.uniform(np.min(sigmas),np.max(sigmas)*3)/np.sqrt(S)
                    elif 'siggam' in mode:
                        mean=-kwargs.get('invader_mean', m.get_param('community_mean') )*(S/Salive)
                        ginv=np.random.choice(gammas)
                        std=np.random.choice([np.min(sigmas),np.max(sigmas)*3])/np.sqrt(Salive)

                    else:
                        mean=-kwargs.get('invader_mean', m.get_param('community_mean') )
                        mstdrel= m.get_param('community_stdrel')
                        if mstdrel:
                            mstd=mstdrel*mean
                        else:
                            mstd= m.get_param('community_std')
                        std=kwargs.get('invader_std',mstd)
                        ginv=kwargs.get('invader_symmetry', m.get_param('community_symmetry') )
                    # print mean,std
                    A0=np.random.normal(0,std,(2,Salive))
                    A0T=A0[::-1,:]
                    if ginv==0:
                        asymm=0
                    else:
                        asymm=(1- np.sqrt(1-ginv**2) )/ginv
                    A0= (A0+asymm * A0T) / np.sqrt(1+ asymm**2 )

                    A0i,Ai0=A0+mean
                else:
                    raise Exception('Use extinct species as invaders: not done yet!')

                if - np.dot(A0i, Nalive)>1 and not 'nonlin' in mode:
                    print 'DIVERGENCE EXPECTED'
                    continue

                # Delta.append(val)
                inv+=1
                result={'mu':mu,'sigma':sigma,'gamma':gamma,'minv':mean,'sinv':std,'ginv':ginv}
                # code_debugger()
                if 'evol' in mode:
                    print "Invader", inv
                    kw={}
                    kw.update(deepcopy(kwargs))
                    kw['n_shape'] = (Salive + 1,)
                    Anew = np.zeros(np.array(Al.shape) + 1)
                    Anew[:-1, :-1] = Al
                    Anew[-1, :-1] = A0i
                    Anew[:-1, -1] = Ai0
                    rnew = np.zeros(Salive + 1)
                    rnew[:-1] = ralive
                    rnew[-1]=r0

                    meps = Model(parameters=deepcopy(prm), dynamics=dyn, **kw)
                    meps.data['growth'].matrix[:] = rnew
                    # code_debugger()
                    if 'foodweb' in mode:
                        # A0i[A0i<0]=0
                        # Ai0[Ai0<0]=0
                        newmats = {}
                        for mat in ('fluxes','efficiency','threshold','community'):
                            matnew=np.zeros(Anew.shape)
                            prev=m.data[mat].matrix[np.ix_(alive,alive)]
                            matnew[:-1,:-1]=prev
                            if mat in ['threshold','efficiency']:
                                matnew[-1]=np.random.choice(prev[prev>0],size=matnew.shape[0])#*ifelse(mat=='threshold',10,1)
                                matnew[:,-1]=np.random.choice(prev[prev>0],size=matnew.shape[0])#*ifelse(mat=='threshold',10,1)
                            else:
                                matnew[-1]=0
                                matnew[:,-1]=0

                            if not 'nontrophic' in mode:
                                if mat=='fluxes':
                                    if r0<0:
                                        matnew[-1,:-1] = np.clip(A0i,0,None)/(0.001+np.mean(A0i<0))  #-np.clip(Ai0, None,0)
                                    else:
                                        matnew[:-1,-1] =np.clip(Ai0,0,None)/(0.001+np.mean(Ai0<0))#-np.clip(A0i, None,0)
                                    # code_debugger()

                            # if mat=='fluxes':
                            #     matnew[-1, :-1] = A0i
                            #     matnew[:-1, -1] = Ai0
                            newmats[mat]=matnew
                            meps.data[mat].matrix[:]=matnew
                        Anew[:-1,:-1]=0
                        if 'nontrophic' in mode:
                            meps.data['nontrophic'].matrix[:] = -Anew
                        meps.data['community'].matrix[:] += -Anew

                    else:
                        meps.data['community'].matrix[:] = -Anew


                    #SMALL INVASION
                    oldN = np.concatenate([Nalive, [0]])
                    x0 = oldN.copy()
                    x0[-1] = 0.01#np.max(x0)*2



                    meps.set_init(n=x0)
                    meps.evol(tmax=tmax, tsample=tsample, converge=1, reseed=False)
                    N1=meps.results['n'][-1]
                    # print np.dot(Anew, x0)[-1] - r0, np.mean(A0i), 1. / np.sum(Nalive), np.dot(A0i, x0[:-1])

                    A00=meps.data['selfint'].matrix[-1]
                    if 'nonlin' in mode:
                        invrate=np.dot(meps.get_jac(0, x0), (x0-oldN))
                        Jac=meps.get_jac(0, x0)
                        Acomp = -np.dot(np.diag(1 / x0), Jac)
                        A00 = -meps.get_jac2(0,x0)[-1,-1]/2.
                        # code_debugger()
                        np.fill_diagonal(Acomp, 0)
                        A0ic=Acomp[-1,:-1]
                        Ai0c = Acomp[:-1, -1]
                        if  invrate[-1]<0 or meps.results['n'][-1,-1]<10**-5:
                            print 'CANNOT INVADE'
                            if 'food' in mode and 0:
                                foodshow(meps)
                                # code_debugger()
                                meps.get_dx(0, x0,frdebug=1)

                                continue

                    rel=(N1/(0.01+oldN) - 1)
                    rel[ np.abs(N1-oldN)<10*death ]=0
                    rel=rel[:-1]
                    result['K0']=r0/A00
                    result['A00']=A00
                    result['resisted']=N1[-1]<death
                    result['N0'] = N1[-1]
                    result['DNtot'] = lanorm(N1 - oldN)
                    result['Nmean']=np.mean(Nalive)
                    result['Nmedian']=np.median(Nalive)
                    #
                    # code_debugger()
                    #STRONG INVASION
                    # meps2 = Model(parameters=deepcopy(prm), dynamics=dyn, **kw)
                    # meps2.data['growth'].matrix[:] = rnew
                    # if 'foodweb' in mode:
                    #     for mat in ('fluxes','efficiency','threshold'):
                    #         meps2.data[mat].matrix[:] = newmats[mat]
                    #     meps2.data['nontrophic'].matrix[:] = -Anew
                    # else:
                    #     meps2.data['community'].matrix[:] = -Anew
                    meps2=meps.copy()
                    x0[-1] = np.max(x0)*5
                    meps2.set_init(n=x0)
                    meps2.evol(tmax=tmax, tsample=tsample, converge=1, reseed=False)
                    N2=meps2.results['n'][-1]
                    rel2=(N2/(0.01+oldN) - 1)
                    rel2[ np.abs(N2-oldN)<10*death ]=0
                    rel2=rel2[:-1]
                    result['resisted_strong']=N2[-1]<death
                    result['N0_strong'] = N2[-1]
                    result['DN_strong'] = lanorm(N2 - oldN)

                    for suffix,rr in zip(('','_obs','_strong','_obs_strong'),(rel,rel[observed],rel2,rel2[observed] )):
                        result['relDNmed'+suffix] = np.median( rr )
                        result['relDNnorm'+suffix] = np.mean( rr**2 )**.5

                    for suffix, rr in zip(('', '_obs', '_strong', '_obs_strong'),
                                              (N1[:-1], N1[:-1][observed], N2[:-1], N2[:-1][observed])):
                        result['lost'+suffix]=lost=np.where(rr<death )[0]
                        result['#lost'+suffix]=len(lost)
                        result['%lost'+suffix]= len(lost)*1. / Salive

                        if len(rr)>len(observed):
                            oldie=Nalive
                            As = A0i
                        else:
                            oldie = Nalive[observed]
                            As = A0i[observed]
                        result['weiDN'+suffix] = np.dot(As, (rr - oldie) )
                        result['weiN'+suffix] = np.dot(As, rr )
                        result['relweiDN'+suffix] = result['weiDN'] / result['weiN']
                        result['DN'+suffix] = lanorm(rr - oldie)

                    # if result['#lost_strong']<result['#lost']:
                    # code_debugger()

                if 'nonlin' in mode and 'evol' in mode:
                    A0i, Ai0 = A0ic, Ai0c
                else:
                    A00 = meps.data['selfint'][-1]

                # code_debugger()
                result['community']=PRMIDX*ncomm+com
                # result.update(measure)
                result['S']= S
                result['S*']= Salive
                result['phi']= Salive*1./S
                result['worstAi0_std']=np.std(worstAi0)
                result['worstAi0_mean']=np.mean(worstAi0)
                result['worstA0i_std']=np.std(worstA0i)
                result['worstA0i_mean']=np.mean(worstA0i)
                result['worstgamma0']=np.corrcoef(worstA0i,worstAi0)[0,1]
                result['worstAi0_IPR']=np.sum(worstAi0**4)
                result['worstA0i_IPR']=np.sum(worstA0i**4)

                for suffix in ('','_obs'):
                    if 'obs' in suffix  and not notobs:
                            continue
                    r0alone=r0
                    if 'obs' in suffix:
                        r0alone += np.dot( A0i[notobs], np.dot(Vnotobs,m.data['growth'].matrix[notobs]  ) )
                    result['r0alone'+suffix] = r0alone
                    result['r0true'+suffix] =r0 -np.dot(A0i, Nalive)
                    result['f0true'+suffix] = -A00 +np.dot(A0i,np.dot(dNdK,Ai0))
                    if 'obs' in suffix:
                        result['f0fixed'+suffix] = -A00+ np.dot(A0i[notobs],np.dot(Vnotobs,Ai0[notobs] ))
                    else:
                        result['f0fixed'+suffix]=-A00
                    result['f0alone'+suffix] =-A00

                result['W'] =result['r0true']

                result['U']=-result['f0true']/np.abs(A00)
                result['alpha0']=np.sqrt(1-np.clip(result['U'],None,0))
                result['V']=1./result['U']
                result['max_AAT']=np.max(A0i*Ai0)
                result['vtrace']=vtrace
                result['vsum']=vsum
                result['Umean']=1-result['minv']**2 * result['vsum'] - result['sinv']**2*result['ginv']*result['vtrace']
                result['sumtraceratio']=np.abs(result['minv']**2 * result['vsum']) / np.abs( result['sinv']**2*result['ginv']*result['vtrace'])
                result['worst']=worst
                try:
                    result['worstinv']=1-lanorm(A0i)* worst * lanorm(Ai0)
                except:
                    code_debugger()
                    continue
                result['prmidx']=PRMIDX

                # code_debugger()
                table.append(result)


    table=pd.DataFrame(table)
    table.to_csv(path+fname+'.csv',index=0)
    return table


def show_worst_mean(table):
    axes=['mu','sigma','gamma','U','V']
    mus,sigmas,gammas,us,vs = [table[a] for a in axes ]
    plt.figure()
    #Worst case
    worsts=table['worst']
    plt.subplot(121)
    plot(xs,np.mean(worsts,axis=1),hold=1,title=r'Worst case $||B^{-1}||$' )

    #Mean case
    us,vs=table['vtrace'],table['vsum']
    plt.subplot(122)
    print np.mean(vs,axis=1)
    print np.mean(us,axis=1)
    plot(xs,mus**2 * np.mean(us,axis=1) + sigmas**2 * gammas * np.mean(vs,axis=1),hold=1 ,
        title=r'Mean case $\mu^2 \sum_{ij} B^{-1}_{ij} + \sigma^2 \gamma \sum_i B^{-1}_{ii} $')


def show_evol(table,**kwargs):
    # code_debugger()

    for suffix in ('','_strong'):
        if not '%lost'+suffix in table:
            table['%lost'+suffix]=table['#lost'+suffix]/table['S*']

    def mkregime(table):
        regime = {}
        regime['no invasion'] = ESS = np.logical_and(table['resisted'] > 0, table['#lost_strong'] == 0)
        regime['coexistence'] = np.logical_and(table['#lost'] == 0, table['#lost_strong'] == 0) & (~ESS)
        if 'nonlin' in kwargs.get('mode', ''):
            regime['coexistence'] = (table['%lost_strong'] < .01) & (~ESS)
        regime['turnover'] = (~regime['coexistence']) & (table['#lost'] > 0) & (table['U'] >= 0)
        regime['irreversible'] = (table['#lost'] > 0) & (table['f0true'] > 0)
        regime['alternative'] = np.logical_and(table['#lost'] == 0, table['#lost_strong'] > 0)
        return regime

    regimes = ['coexistence', 'no invasion', 'turnover', 'irreversible', 'alternative']
    colors = ['b', 'r', 'k', 'g', 'y']


    truetable = table
    for i, table in truetable.groupby('prmidx'):
        plt.figure()
        plt.suptitle('Mu {} sigma {}'.format(table['mu'].median(), table['sigma'].median()))
        # plt.subplot(131)
        # scatter(us,vs,xlabel=r'Tr($B^{-1}$)', ylabel=r'$\sum_{ij} B^{-1}_{ij}$',c=table["#lost"]*1./table['S*'],hold=1)# ,cmap=plt.cm.Reds)
        for axes in [('f0true', 'r0true')]:  # ('alpha0','K0')]:
            plt.subplot(121)
            regime = mkregime(table)
            legends = []
            ngrid = 300
            xs, ys = table[axes[0]], table[axes[1]]
            if 'U' in axes :
                xmin, xmax = -1.25, 1.25
                ymin, ymax = -1.25, 1.25
            else:
                xmin, xmax = np.min(xs), np.max(xs)
                ymin, ymax = np.min(ys), np.max(ys)
            X, Y = np.meshgrid(np.linspace(xmin, xmax, ngrid),
                               np.linspace(ymin, ymax, ngrid))
            XY = np.vstack([X.ravel(), Y.ravel()]).T
            ZS = []
            # MAT=np.zeros(X.shape+(3,) )
            idx = -1
            xtmp,ytmp,regtmp=[],[],[]
            for ir,i in enumerate(regimes):
                j=regime[i]
                idx += 1
                if np.sum(j) == 0:
                    ZS.append(None)
                    continue
                legends.append(i)
                xs, ys = table[axes[0]][j], table[axes[1]][j]
                xs, ys = xs.values, ys.values
                nx=list(xs)[:300]
                xtmp+=nx
                ytmp+=list(ys)[:len(nx)]
                regtmp+=[ir for nn in nx ]
                from sklearn.neighbors.kde import KernelDensity
                def bandwidth(X):
                    return 1.06 * np.std(X) / (X.shape[1]) ** 0.2

                try:
                    mat = np.array((np.array(xs), np.array(ys)))
                except:
                    code_debugger()
                # kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth(mat)).fit(mat.T)
                # Z = np.exp(kde.score_samples(XY)).reshape(X.shape)
                # ZS.append(Z)
                # colors= [np.array((1.,0,0)),np.array((0,1.,0)),np.array((0,0,1.)),np.array((1.,1.,1.)) ]
                # MAT+= np.multiply.outer(Z/np.max(Z),colors[idx])  /len(regime)
                scatter(xs[:1], ys[:1], xlabel=axes[0], ylabel=axes[1], c=colors[idx], hold=1, newfig=0, s=15, alpha=0.5)
            plt.legend(legends)
            regs = np.argmax(np.array([regime[j].values for j in regimes]),axis=0)
            from sklearn import svm
            xs, ys = table[axes[0]], table[axes[1]]
            # select=(xs>xmin)&(xs<xmax)&(ys>ymin)&(ys<ymax)
            # xs,ys,regs=xs[select],ys[select],regs[select]
            # code_debugger()
            xs,ys,regs=xtmp,ytmp,regtmp

            try:
                from sklearn.model_selection import GridSearchCV
                parameters = {'kernel': ['linear', 'rbf'], 'C': [.1, 1., 10.], 'gamma': [.2, .5, 0.7, 1.]}
                svc = svm.SVC(gamma=0.7, kernel='rbf', C=1.)
                clf = GridSearchCV(svc, parameters, cv=5)
                SVM = clf.fit(np.array([np.array(xs), np.array(ys)]).T, regs)
                Z = SVM.predict(XY).reshape(X.shape)
            except Exception as e:
                print e
                Z = np.ones(X.shape)

            # code_debugger()
            # plt.contourf(X, Y, Z, colors=colors, alpha=.3,levels=np.array(range(len(np.unique(Z))+1))-.5 )


            for idx,i in enumerate(regimes):
                plt.contourf(X, Y, (Z==idx).astype('int'), colors=[colors[idx]], alpha=.25,levels=(.5,1.5))#np.array(range(len(np.unique(Z))+1))-.5 )

            "ONLY PLOT FIRST 200 POINTS IN EACH REGIME"
            for idx,i in enumerate(regimes):
                j=regime[i]
                xs, ys = table[axes[0]][j], table[axes[1]][j]
                xs, ys = xs.values, ys.values
                scatter(xs[:200], ys[:200], xlabel=axes[0], ylabel=axes[1], c=colors[idx], hold=1, newfig=0, s=15, alpha=0.6)


            if 'f0true' in axes:
                plt.axhline(0, color='k')
                plt.axvline(0, color='k')
            # plt.subplot(122)
            # plt.imshow(MAT[::-1]/np.max(MAT))
            # plt.axhline(0, color='w')
            # plt.axvline(0, color='w')
            plt.xlim(xmin=xmin, xmax=xmax),plt.ylim(ymin=ymin,ymax=ymax)

            # plt.gca().set_aspect('equal')
            plt.subplot(122)
            xs, ys = table[axes[0]].values, table[axes[1]].values
            select=[j for i in regimes[::-1] for j in np.where(regime[i])[0][:200]]
            handle = scatter(xs[select], ys[select], c=(table['%lost_strong'] ).values[select], cmap=plt.cm.Reds, hold=1,
                             xlabel=axes[0],s=5, ylabel=axes[1], title='Fraction lost')  # ,plt.colorbar()
            # plt.gca().set_aspect('equal')
            plt.colorbar(handle, ax=plt.gca())
            if 'f0true' in axes:
                plt.axhline(0, color='k')
                plt.axvline(0, color='k')
            plt.xlim(xmin=xmin, xmax=xmax), plt.ylim(ymin=ymin, ymax=ymax)

    regime=mkregime(truetable)

    init=['','_strong'][0]
    scaling=''#['','rel','N0','relN0'][0]
    Nrel =  np.abs(table['W']*table['V'])/table['K0']
    Nrel[regime['irreversible']]=1


    # if 'N0' in scaling:
    #     Nrel = table['N0' + init] / table.K0
    # table['relDNmed2' + init] = table['relDNmed' + init]
    # table['relDNnorm2' + init] = table['relDNnorm' + init]
    # impact2 = (1. / table['V'] - 1) / (table['W'] / table['K0'] - 1) * Nrel
    # impact2=(table['U']-1)/(table['W']-table['K0'])*np.abs(table['V']*table['W'])
    # if 'rel' in scaling:
    #     impact2 /= Nrel
    #     table['relDNmed2' + init] /= Nrel
    #     table['relDNnorm2' + init] /= Nrel
    for suffix in  ('','_obs','_strong'):
        table=truetable
        if not 'f0fixed'+suffix in table:
            continue
        N0ref=np.abs( table['r0true']/table['f0true'] )
        # N0ref=table['N0' + init]
        N0ref[regime['irreversible']]=table['K0'].values[regime['irreversible']]
        impact=np.abs( (table['f0true']-table['f0fixed'+suffix] )/(table['r0true']-table['r0alone'+suffix] ) ) *N0ref
        table['impact'+suffix]=impact
        plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize']) * (1.2, 1.6))
        plt.suptitle(suffix)
        if not '%lost'+suffix in table:
            table['%lost' + suffix]=table['%lost']
        xs,ys,zs,ex=np.log10(np.abs(-table['impact'])),np.log10(np.abs(-table['relDNmed'+suffix])),np.log10(table['relDNnorm'+suffix]),table['%lost'+suffix]
        good = (~np.isnan(xs)) & (~np.isnan(ys)) & (~np.isnan(zs))&(~np.isinf(xs))&(~np.isinf(ys))&(~np.isinf(zs)) &(zs<4)
        xs,ys,zs,ex=xs[good],ys[good],zs[good],ex[good]
        rg = min(np.min(xs), np.min(ys)), max(np.max(xs), np.max(ys))
        rg2 = min(np.min(xs), np.min(zs)), max(np.max(xs), np.max(zs))
        from scipy.stats import spearmanr
        rex=spearmanr(xs,ex)[0]
        rmed=spearmanr(xs,ys)[0]
        rtot=spearmanr(xs,zs)[0]
        # print xs,ys,zs
        if 0:
            plt.subplot(144)
            bins=np.linspace(-1,len(xs)+1,3)
            xbins=np.digitize(np.argsort(np.argsort(xs)),bins)
            nbins=set(xbins)
            xmean=[np.mean(xs[xbins==b]) for b in nbins]
            for ws in (ex,ys,zs):
                tmp=[spearmanr(xs[xbins==b], ws[xbins==b])[0] for b in nbins]
                plt.plot(xmean,tmp)
            plt.legend(['Extinctions','Median','Absolute'])
        plt.subplot(222)
        plot(rg, rg, hold=1,color='k')
        plt.subplot(223)
        plot(rg2,rg2,hold=1,color='k')
        # plt.subplot(224)
        # plot(rg2,rg2,hold=1,color='k')
        for col,reg in zip(colors,regimes):
            if reg=='no invasion' or (init=='' and reg=='alternative'):
                continue
            table=truetable[regime[reg]]

            plt.subplot(221)
            plt.xlim(xmin=-5,xmax=3)
            plt.ylim(ymin=-0.1,ymax=1.1)
            xs=np.log10(np.abs(-table['impact'+suffix]))
            try:
                scatter(xs,table['%lost'+suffix],alpha=.5,c=col,log='',hold=1,xlabel='log10 Predicted impact',
                    ylabel='Extinction fraction',title='Extinctions, rho={:.2g}'.format(rex))
            except:
                pass
            plt.subplot(222)
            ys=np.log10(np.abs(table['DN'+suffix])/(table['Nmean']))
            ys=np.log10(np.abs(table['relweiDN'+suffix]))
            # print spearmanr(xs,ys)
            scatter(xs,ys,log='',alpha=.5,c=col,xlabel=r'$\log_{10}$ Predicted impact',ylabel='log10 Weighted impact',
                    title='Median impact, rho={:.2g}'.format(rmed) ,hold=1)
            plt.subplot(223)
            xs,zs=np.log10(np.abs(-table['impact'+suffix])),np.log10(table['relDNnorm'+suffix])
            # xs,zs=xs[zs<4],zs[zs<4]
            plt.xlim(xmin=-5,xmax=3)
            plt.ylim(ymin=-5,ymax=3)
            scatter(xs,zs,alpha=.5,c=col,xlabel=r'$\log_{10}$ Predicted impact',ylabel=r'$\log_{10} \;|| \Delta N/N ||$',
                    hold=1,title='Abundance change, rho={:.2g}'.format(rtot ))
            plt.subplot(224)
            xs,zs=np.log10(table['r0true']),np.log10(np.abs(table['relweiDN'+suffix]))
            plt.ylim(ymin=-5,ymax=3)
            rtot=spearmanr(xs[table['r0true']>0],zs[table['r0true']>0])[0]
            scatter(xs,zs,alpha=.5,c=col,xlabel=r'$\log_{10} r_0^{\mathrm{true}}$',ylabel=r'$\log_{10} \;|| \Delta N/N ||$',
                    hold=1,title='Abundance change, rho={:.2g}'.format(rtot ))
    plt.show()
    code_debugger()


def notassembled(ntrials=500,S=80):
    import seaborn as sns
    res = 50
    sigmas = np.linspace(0.05, .8, res)
    gammas = np.array([-.99, 0, .99])
    sigmas, gammas = np.multiply.outer(sigmas, np.ones(gammas.shape)), np.multiply.outer(np.ones(sigmas.shape), gammas)
    sigmas, gammas = sigmas.ravel(), gammas.ravel()
    mus = np.ones(sigmas.shape) * 5.

    table=[]
    for i in range(ntrials):
        mn,sd,g=np.random.choice(mus)/S,np.random.choice(sigmas)/np.sqrt(S),np.random.choice(gammas)
        M= np.random.normal(0,sd, (S,S))
        np.fill_diagonal(M,0)
        if g == 0:
            asymm = 0
        else:
            asymm = (1 - np.sqrt(1 - g ** 2)) / g
        M = (M + asymm * M.T) / np.sqrt(1 + asymm ** 2) + mn
        # print g, np.corrcoef(offdiag(M),offdiag(M.T))[0,1]

        I = np.eye(S)
        from scipy.linalg import inv as lainv, norm as lanorm
        A=I - M


        dNdK = lainv(A)
        worst = lanorm(dNdK, ord=2)

        result={'sigma':sd*np.sqrt(S),'gamma':g,'worst': worst}
        result['vtrace'] = np.trace(dNdK)
        result['vsum'] = np.sum(dNdK)

        table.append( result )
    table=pd.DataFrame(table)
    sns.relplot(x='sigma',y='worst',col='gamma',data=table)
    sns.relplot(x='sigma',y='vtrace',col='gamma',data=table)
    sns.relplot(x='sigma',y='vsum',col='gamma',data=table)
    plt.show()


def smooth(X,Y,Z,res=20):
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.model_selection import GridSearchCV
    clf = GridSearchCV(KernelRidge(kernel='rbf', gamma=0.1), cv=5,
                       param_grid={"alpha": [1e0, 0.1, 1e-2, 1e-3],
                                   "gamma": np.logspace(-2, 2, 5)})
    XYsmooth = np.array([x for x in itertools.product(np.linspace(np.min(X), np.max(X), res),
                                                      np.linspace(np.min(Y), np.max(Y), res),)])
    #
    #
    XY=np.array([X,Y]).T
    clf.fit(XY, Z)
    print 'Fitted'
    Zsmooth = clf.predict(XYsmooth)
    return XYsmooth,Zsmooth


def show_sigmadiff(table,hold=0):
    import seaborn as sns
    medS=int(np.round(table["S"].median()))
    table['vtrace']/=table['S*']
    table['vsum']/=table['S*']**2
    table['sumtraceratio']=np.log10(table['sumtraceratio'])

    hue=['sinv','mu'][-1]



    table['sg']=table['sigma']**2 * table['gamma']
    table['sgp']=table['sg'] *table['phi']
    table['sg0'] = table['sinv'] ** 2 * table['ginv']
    table['worstheo']=1/( 1-table['mu']/table['S'] -(table['gamma']+1)*table['phi']*table['sigma']**2 )
    table['worstinvtheo']=1-(1-table['worstinv'])/table['worst'] * table.worstheo
    scatter(table.worstheo,table.worst)

    if 1:
        tab=table[table.mu<2].copy()
        fig, axes = plt.subplots(ncols=1, figsize=np.array(mpfig.rcParams['figure.figsize']) * (1., 1))
        ax=axes
        ax.fill_between([-1,1],-100,color='k',alpha=.1)
        sns.scatterplot(x='sgp',y='U',data=tab.iloc[np.random.random(tab.shape[0] )<.01],hue='sg0',ax=ax )
        sg0bins=np.linspace(tab.sg0.min(),tab.sg0.max(),4)
        plt.legend([r'$\sigma^2_0 \gamma_0=${:.2g}'.format(f) for f in np.mean([sg0bins[1:],sg0bins[:-1]],axis=0)  ])
        for sg0,gp in tab.groupby(np.digitize(tab.sg0,sg0bins)):
            gp=gp.groupby('sgp').mean()
            # plt.plot(gp.index,gp.Umean,linestyle='--',lw=2,color='g' )
            if len(gp.index)<20:
                continue
            from scipy.signal import savgol_filter
            ys=savgol_filter(gp.Umean.values,11,3)
            plt.plot(gp.index, ys, linestyle='--', lw=3, color='k')
        # plt.fill_between( )
        gp=tab[tab.sg0<-.05].groupby('sg').min()
        # plt.plot(gp.index,gp.worstinv,color='k',lw=3 )
        plt.plot(gp.index,gp.worstinvtheo,color='k',lw=3 )
        plt.xlabel(r'$\sigma^2 \gamma$')
        plt.show()
    code_debugger()
    sns.relplot(x='sg',y='sg0',data=table,hue='vtrace' )
    sns.relplot(x='sg',y='sg0',data=table,hue='vsum' )
    sns.relplot(x='sg',y='sg0',data=table,hue='worst' )
    sns.relplot(x='sg',y='sg0',data=table,hue='U' )



    if 1:
        plt.figure()
        plt.colorbar(plt.scatter(table['sg'],table['sg0'], c=table['U'],cmap='seismic_r',vmax=5,vmin=-5))
        plt.show()

    sns.relplot(x='sigma',y='worst',col='gamma',data=table,hue=hue)

    sns.relplot(x='sigma',y='worst',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='vtrace',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='vsum',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='sumtraceratio',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='S*',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='worstgamma0',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='worstAi0_mean',col='gamma',data=table,hue=hue)
    #table['worstAi0_std_rescaled']=table['worstAi0_std']*np.sqrt(table['S*'])
    #sns.relplot(x='sigma',y='worstAi0_std',col='gamma',data=table,hue=hue)
    #sns.relplot(x='sigma',y='worstAi0_std_rescaled',col='gamma',data=table,hue=hue)
    table['worstAi0_IPR_rescaled']=table['worstAi0_IPR']*(table['S*'])
    sns.relplot(x='sigma',y='worstAi0_IPR',col='gamma',data=table,hue=hue)
    sns.relplot(x='sigma',y='worstAi0_IPR_rescaled',col='gamma',data=table,hue=hue)
    if 0:
        sns.relplot(x='sigma',y='sinv',row='gamma',col='ginv',data=table[table['mu']==1],hue='Umean')
        sns.relplot(x='sigma',y='sinv',row='gamma',col='ginv',data=table[table['mu']==1],hue='worst')
        sns.relplot(x='sinv',y='Umean',row='gamma',col='ginv',data=table,hue='sigma')
        sns.relplot(x='sinv',y='Umean',row='gamma',col='ginv',data=table,hue='mu')

    plt.figure()
    mu1=table['mu'].min()
    plt.subplot(221)
    vmin,vmax=table['Umean'].min(),table['Umean'].max()
    tab=table[(table['mu']==mu1) & (table['gamma']>.5)&(table['ginv']>.5) ]
    plt.scatter(tab['sigma'],tab['sinv'],c=tab['Umean'].values,cmap='seismic_r',vmin=vmin,vmax=vmax),plt.colorbar(),plt.xlabel('sigma'),plt.ylabel('sinv'),plt.title('Umean (ginv=1)')
    plt.subplot(222)
    tab=table[(table['mu']==mu1) & (table['gamma']>.5)&(table['ginv']<-.5) ]
    plt.scatter(tab['sigma'],tab['sinv'],c=tab['Umean'].values,cmap='seismic_r',vmin=vmin,vmax=vmax),plt.colorbar(),plt.xlabel('sigma'),plt.ylabel('sinv'),plt.title('Umean (ginv=-1)')

    plt.subplot(223)
    table['logworst']=np.log10(-table['worstinv'])
    tab=table[(table['mu']==mu1) & (table['gamma']>.5)&(table['ginv']>.5) ]
    vmin,vmax=tab['logworst'].min(),tab['logworst'].max()
    plt.scatter(tab['sigma'],tab['sinv'],c=tab['logworst'].values,vmin=vmin,vmax=vmax,cmap='Reds'),plt.colorbar(),plt.xlabel('sigma'),plt.ylabel('sinv'),plt.title('worst (gamma=1)')
    plt.subplot(224)
    tab=table[(table['mu']==mu1) & (table['gamma']<-.5)&(table['ginv']>.5) ]
    plt.scatter(tab['sigma'],tab['sinv'],c=tab['logworst'].values,vmin=vmin,vmax=vmax,cmap='Reds'),plt.colorbar(),plt.xlabel('sigma'),plt.ylabel('sinv'),plt.title('worst (gamma=-1)')
    #plt.show()
    # code_debugger()


    table['cUmean']=1-table['Umean']
    table['cU']=1-table['U']
    tab=table.groupby(['sigma','gamma','ginv']).median().reset_index()
    # sns.relplot(x='worst',y='cU',row='gamma',col='ginv',data=tab,hue='sigma')
    sns.relplot(x='cUmean',y='cU',row='gamma',col='ginv',data=table,hue=hue)
    # scatter(table['max_AAT'],table['cU'],newfig=1,hold=1,xlabel='max(A0i*Ai0)',title='1-U')
    sns.relplot(x='sumtraceratio',y='cU',row='gamma',col='ginv',data=table,hue=hue)
    if not hold:
        plt.show()


def foodweb(S=32,ncomm=30,ninv=30):
    prm = {
        'n': {
            'type': 'variable',
            'axes': ('com',),
            'death': 10. ** -9,
            'shape': (S,),
            'mean': 1.,
            'std': 1,
            'sign': 1,
        },
        'growth': {
            'type': 'matrix',
            'variables': ['n'],
            'requires': ['niches', ],
            'dynamics': 'nlv',
            'role': 'growth',
            'structure': 'niche',

        },

        'nontrophic': {
            'type': 'matrix',
            'role': 'nontrophic',
            'dynamics': 'nlv',
            'variables': [('n', 'com'), ('n', 'com')],
            'mean': 0,
            'std':0,
        },
        'selfint': {
            'type': 'matrix',
            'variables': ['n'],
            'dynamics': 'nlv',
            'role': 'diagonal',
            'mean': 1, 'std': 0.,
            'sign': 1,
            'structure': 'trophic',
            'requires': ['niches'],
        },
        'fluxes': {
            'type': 'matrix',
            'variables': [('n', 'com'), ('n', 'com')],
            'role': 'fluxes',
            'requires': ['niches', 'edges'],
            'structure': 'niche',
            'dynamics': 'nlv',
            'diagonal': 0,
            'sign': 1,
            'edgefrom': 'edges',
        },
        'community': {
            'type': 'matrix',
            'variables': [('n', 'com'), ('n', 'com')],
            'role': 'interactions',
            'structure': 'trophic',
            'requires': ['fluxes', 'efficiency', 'niches'],
            'dynamics': 'nlv',
            'ordered': 0,
            'diagonal': 0,
            'mean':1,
        },
        'niches': {
            'type': 'constant',
            'role': 'niche',
            'axes': [('n', 'com')],
            'distribution': 'randint',
            'min': 1,
            'max': 3,
            'save': True,
        },
        'edges': {
            'type': 'matrix',
            'variables': [('n', 'com'), ('n', 'com')],
            'mode': 'network',
            'connectivity': .4,
            'role': 'edges',
            'save': False,
            'structure': 'random',
        },

        'efficiency': {
            'type': 'matrix',
            'role': 'efficiency',
            'variables': ['n', 'n'],
            'dynamics': 'nlv',
            'structure': 'niche',
            'requires': ['niches'],
            'mean': 1,
            'std': 0,
            'sign': 1,
            'factor': 1,
        },

        'threshold': {
            'type': 'matrix',
            'variables': [('n', 'com'), ('n', 'com')],
            'dynamics': 'nlv',
            'role': 'threshold',
            'mean': .2,
            'stdrel': 0,
        },

        'exponent': {
            'type': 'matrix',
            'variables': ['n'],
            'dynamics': 'nlv',
            'role': 'exponent',
            'mean': 1,
            'stdrel': 0,
        },
    }

    dyn = {
        'nlv': {'type': [ 'trophicfr'][0],
                'variables': [('n', 'com'), ],
                },
    }

    mode='evol_rndinv_foodweb_obs' #+ '_nonlin'
    df = invasibility(path='criterion_foodweb' , extend='extend' in sys.argv,rerun='rerun' in sys.argv,
                      ncomm=ncomm, ninv=ninv,
                      mode=mode,S=S,prm=prm,dyn=dyn)
    show_evol(df,mode=mode)

if __name__ =='__main__':
    if True in [ 'vw' in x for x in sys.argv]:
        df=invasibility(path='criterion_VW'+ifelse('vw2' in sys.argv,'2',''),rerun= 'rerun' in sys.argv,ncomm=50,
                        mode='evol_rndinv',ninv=50)
        # code_debugger()
        if not 'S*' in df:
            df['S*']=df['S']
            df['S']=df['S']+df['#lost']
        show_evol(df)

    if True in [ 'obs' in x for x in sys.argv]:
        df=invasibility(path='criterion_obs',rerun= 'rerun' in sys.argv,ncomm=40,
                        mode='evol_rndinv_obs',ninv=40,S=100)
        show_evol(df)

    if 'sig' in sys.argv:
        df=invasibility(path='criterion_sigma'+ifelse('2' in sys.argv,'2',''),rerun= 'rerun' in sys.argv,ncomm=1,mode='sigmadiff',ninv=30)
        show_sigmadiff(df)
    if 'sigg' in sys.argv:
        df=invasibility(path='criterion_siggam'+ifelse('2' in sys.argv,'2',''),rerun= 'rerun' in sys.argv,ncomm=30,mode='siggam',ninv=100)
        show_sigmadiff(df)

    if 'food' in sys.argv :
        foodweb(S=100,ncomm=32,ninv=32)

