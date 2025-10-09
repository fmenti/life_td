import numpy as np #arrays
import astropy as ap #votables
import importlib #reloading external functions after modification
import seaborn as sns
import sys  
from matplotlib import pyplot as plt

sys.path.append('../../implementation/life/data_generation/life_td_data_generation') 

#self created modules
from utils.io import save, load, stringtoobject, Path
from utils.catalog_comparison import create_common
import analysis as la
importlib.reload(la)#reload module after changing it
import LIFE_StarCat4 as ltc4
importlib.reload(ltc4)#reload module after changing it
import provider as p
importlib.reload(p)#reload module after changing it

para3=['sim_ra', 'sim_dec', 'sim_plx', 'distance', 'gal_coord_l', 'gal_coord_b', 
       'mod_Teff', 'mod_R', 'mod_M', 'sep_phys',  'sim_i', 'sim_j']
para4=['coo_ra', 'coo_dec', 'plx_value', 'dist_st_value', 'coo_gal_l', 'coo_gal_b', 
       'teff_st_value', 'radius_st_value','mass_st_value','sep_phys_value',  'mag_i_value', 'mag_j_value']

para_hwo=['ra', 'dec', 'sy_plx', 'sy_dist', 'st_teff', 'st_rad', 'st_mass']
para_starcat4=['coo_ra', 'coo_dec', 'plx_value', 'dist_st_value', 'teff_st_value',
               'radius_st_value','mass_st_value']


def analysis(StarCat4,distance_cut_in_pc):
    """
    This Function analyses the StarCat4.
    
    :param StarCat4:
    :type StarCat4:
    :param float distance_cut_in_pc:
    """
    
    la.spechistplot([StarCat4['class_temp'][np.where(StarCat4['binary_flag']=='False')], \
              StarCat4['class_temp'][np.where(StarCat4['binary_flag']=='True')],StarCat4['class_temp']]\
             ,['Singles','Binaries','Total'],'spechist')

    la.final_plot([StarCat4['class_temp','dist_st_value']],['LIFE-StarCat4'],distance_cut_in_pc)

    within45deg=StarCat4[np.where(StarCat4['ecliptic_pm45deg']=='True')]
    
    print('if we cut for architecture reason everything above and below +45 deg, -45 deg respectively')
    print(f'we are left with {np.round(100*len(within45deg)/len(StarCat4),1)}% of our stars')
    print(f'this corresponds to a factor of {np.round(len(StarCat4)/len(within45deg),1)} fewer')
    remaining_singles=len(within45deg[np.where(within45deg['binary_flag']=='False')])/len(within45deg)
    print('from the remaining stars almost all(',np.round(100*remaining_singles,1),'%) are singles')

    return

def create_updated_ltc3():
    """
    This function updates the StarCat3 with newer main_id's.
    
    :returns: LIFE-StarCat3 with updated main_id from SIMBAD.
    :rtype: astropy.table.table.Table
    """
    
    #load LTC3
    ltc3=ap.io.ascii.read("data/LIFE-StarCat3.csv")
    ltc3=stringtoobject(ltc3,3000) #len 1732

    ltc3_update=p.fetch_main_id(ltc3,'sim_name',name='main_id',oid=False)
    #print(ltc3_update)#1728 -> lost 4 objects
    lost=ltc3[np.where(np.invert(np.in1d(ltc3['sim_name'],ltc3_update['sim_name'])))]
    #L  326-61, L  380-78,  L  755-19, LP  399-68 -> maybe were truncated?

    print('To compare the version 3 to version 4 we had to updated the simbad main identifiers of v3.')
    print('Of the originally ',len(ltc3), 'objects ',len(lost),'could not be matched that way.')
    print('After matching on position those could be retrieved again.')
    newltc3=ltc3
    newltc3['sim_name'][np.where(newltc3['sim_name']=='L  326-61')]='L 326-21'
    newltc3['sim_name'][np.where(newltc3['sim_name']=='L  380-78')]='L 380-8'
    newltc3['sim_name'][np.where(newltc3['sim_name']=='L  755-19')]='LP 755-19'
    newltc3['sim_name'][np.where(newltc3['sim_name']=='LP  399-68')]='L 399-68'
    newltc3=p.fetch_main_id(newltc3,'sim_name',name='main_id',oid=False)
    print('newltc3',newltc3)
    save([newltc3],['newltc3'],location='data/')
    [StarCat4]=load(['StarCat4'],location='data/')
    in_3_but_not_in_4=newltc3[np.where(np.invert(np.in1d(newltc3['main_id'],StarCat4['main_id'])))]
    print('in_3_but_not_in_4',in_3_but_not_in_4)
    save([in_3_but_not_in_4],['in_3_but_not_in_4'],location='data/')    
    return newltc3


def prepare_dataframe_for_sns(table,paras,labels):
    """
    This function prepares the dataframe for seaborn boxplot.
    
    :param table: Table containing two data sets and the column 
        names given in paras as well as class_temp column.
    :type table: astropy.table.table.Table
    :param paras: Column names of table with first list corresponding to first lables entry.
    :type paras: list(list(str))
    :param labels: Touple with lable for the two data sets.
    :type labels: list(str)
    :returns: Dataframe ready to be ingested by seaborn boxplot.
    :rtype: pandas.core.frame.DataFrame
    """
    seaborn_join=table[['class_temp']+paras[0]].copy()
    seaborn_join['catalog']=[labels[0] for j in range(len(seaborn_join))]
    temp=table[['class_temp']+paras[1]].copy()
    temp['catalog']=[labels[1] for j in range(len(temp))]
    seaborn_join.rename_columns(paras[0],paras[1])
    seaborn_join=ap.table.vstack([seaborn_join,temp])
    df=seaborn_join.to_pandas()
    return df

def snsplot(df,y,path):
    """
    This function creates a boxplot.
    
    :param df:
    :type df:
    :param y:
    :type y:
    :param str path:
    """
    
    plt.figure()
    sns.boxplot(data=df,x="class_temp",y=y ,hue='catalog',order=['F','G','K','M'])
    plt.savefig(path+'.png')
    plt.show()
    plt.figure(figsize=(2,3))
    sns.boxplot(data=df,y=y ,hue='catalog',legend=False)
    plt.savefig(path+'total.png', bbox_inches = 'tight')
    plt.show()
    return

def threecatboxplot(data,para):
    """
    This function creates a boxplot.
    
    :param data:
    :type data:
    :param para:
    :type para:
    """
    
    fig = plt.figure()
 
    # Creating axes instance
    ax = fig.add_axes([0, 0, 1,1])
    plt.ylabel(para)
    # Creating plot
    bp = ax.boxplot(data, labels=['StarCat3','StarCat4','SPORES'])
    plt.savefig('plots/'+para+'allcat_total.png', bbox_inches = 'tight')
    # show plot
    plt.show()
    return


def create_hwo():
    """
    Prepares hwo catalog for comparison to LIFE-StarCat4.
    
    Reads the HWO SPORES catalog file, adds SIMBAD main identifiers and transforms
    StarCat4 comparable columns from type string to float.
    
    :returns: Table containing HWO catalog data.
    :rtype: astropy.table.table.Table
    """
    
    spores_raw=ap.io.ascii.read("data/spores_catalog_v1.2.0.csv")#,header_start=0,)
    spores_raw=stringtoobject(spores_raw,3000)
    colnames=list(spores_raw[0])
    colnames=colnames[1:]
    data=spores_raw[1:].copy()

    data.remove_columns('col1')
    spores=ap.table.Table(data,names=colnames)
    print(spores)
    #add main_id column
    print(spores.colnames)
    temp=spores['tic_id','hip_name'].copy()
    hwo=p.fetch_main_id(temp,'tic_id',name='main_id',oid=False)
    hwo.rename_column('hip_name','temp')
    hwo=ap.table.join(spores,hwo,keys='tic_id')
    hwo.remove_column('temp')
    print(hwo)
    for para in para_hwo:
            hwo[para]=hwo[para].astype(float)
    save([hwo],['hwo'],location='data/')
    return hwo

def create_hpic():
    HPIC= ap.io.ascii.read('data/full_HPIC.csv')
    print(HPIC)
    print(HPIC.colnames)
    HPIC=HPIC[np.where(HPIC['sy_dist']!='null')]#143 null values in distance
    HPIC=HPIC[np.where(HPIC['st_spectype']!='null')]#of those about 2600 have no spectral type given
    HPIC['sy_dist']=HPIC['sy_dist'].astype(float)
    save([HPIC],['hpic'],location='data/')
    return HPIC


def distance_boxplot_catalog_comparison(ax2):
    
    #don't have a join thingy so might need to adapt stuff
    #df=prepare_dataframe_for_sns(join_34,[para3,para4],labels)
    #I need class_temp for spectral type sorting. can do that too with the specsomething function
    
    seaborn_join=table[['class_temp']+paras[0]].copy()
    seaborn_join['catalog']=[labels[0] for j in range(len(seaborn_join))]
    temp=table[['class_temp']+paras[1]].copy()
    temp['catalog']=[labels[1] for j in range(len(temp))]
    seaborn_join.rename_columns(paras[0],paras[1])
    seaborn_join=ap.table.vstack([seaborn_join,temp])
    df=seaborn_join.to_pandas()
    
    snsplot(df,'dist_st_value')
    
    
    return 

def ltc4_comparison(create=False,distance_cut_in_pc=30.,boxplots=True,finalplots=True):
    """
    This function compares StarCat4 to StarCat3 and HWO catalog SPORES.
    
    :param bool load: Defaults to True, determines if tables are created or loaded.
    :param float distance_cut_in_pc: Defaults to 30., determines distance cut of StarCat4 if load=False.
    """
    
    if create:
        StarCat4=ltc4.StarCat_creation(distance_cut_in_pc,querying=True)
        newltc3=create_updated_ltc3()
        common_newltc3,common_StarCat43=create_common(newltc3,StarCat4)
        join_34=ap.table.join(common_newltc3,common_StarCat43,keys='main_id')
        save([join_34],['join_34'],location='data/')   
        hwo=create_hwo()
        common_hwo,common_StarCat4hwo=create_common(hwo,StarCat4)
        save([common_newltc3,common_StarCat43,
                common_hwo,common_StarCat4hwo],
                ['common_newltc3','common_StarCat43',
                'common_hwo','common_StarCat4hwo'],location='data/')   
        join_hwo4=ap.table.join(common_hwo,common_StarCat4hwo,keys='main_id')
        save([join_hwo4],['join_hwo4'],location='data/')  
        
        hpic=create_hpic()
        
    else:
        [StarCat4]=load(['StarCat4'],location='data/')
        [newltc3]=load(['newltc3'],location='data/')
        [hpic]=load(['hpic'],location='data/')
        
        [join_34]=load(['join_34'],location='data/')
        [join_hwo4]=load(['join_hwo4'],location='data/')
        [hwo]=load(['hwo'],location='data/')
        [common_newltc3,common_StarCat43,common_hwo,common_StarCat4hwo]=\
                load(['common_newltc3','common_StarCat43',
                'common_hwo','common_StarCat4hwo'],location='data/')
        
           
    if boxplots:
        #comparing only objects that are in both catalogs
        df_34=prepare_dataframe_for_sns(join_34,[para3,para4],['3','4'])

        for y in para4:
            snsplot(df_34,y,'plots/sns_3_4'+y)

        df_hwo4=prepare_dataframe_for_sns(join_hwo4,[para_hwo,para_starcat4],['SPORES','StarCat4'])

        for y in para_starcat4:
            snsplot(df_hwo4,y,'plots/sns_hwo_4'+y)
    

        #comparing all objects in the catalogs
        data = [newltc3['distance'],StarCat4['dist_st_value'], hwo['sy_dist']]
        threecatboxplot(data,'Stellar_Distance')

        data = [newltc3['mod_M'],StarCat4['mass_st_value'], hwo['st_mass']]
        threecatboxplot(data,'Stellar_Mass')

        data = [newltc3['mod_R'],StarCat4['radius_st_value'], hwo['st_rad']]
        threecatboxplot(data,'Stellar_Radius')

        data = [newltc3['mod_Teff'],StarCat4['teff_st_value'], hwo['st_teff']]
        threecatboxplot(data,'Stellar_Effective_Temperature')
    
    if finalplots:
    
        la.final_plot([StarCat4['class_temp','dist_st_value'],newltc3['sim_sptype','distance']],
                  ['LIFE-StarCat4','LIFE-StarCat3'],path='plots/final_plot_34.png')
    
        la.final_plot([common_StarCat43['class_temp','dist_st_value'],common_newltc3['sim_sptype','distance']],
                  ['common_LIFE-StarCat4','common_LIFE-StarCat3'],path='plots/final_plot_common_34.png')
    
        la.final_plot([StarCat4['class_temp','dist_st_value'],hwo['st_spectype','sy_dist']],
                  ['LIFE-StarCat4','SPORES'],path='plots/final_plot_hwo.png')
    
        la.final_plot([common_StarCat4hwo['class_temp','dist_st_value'],common_hwo['st_spectype','sy_dist']],
                  ['LIFE-StarCat4','SPORES'],path='plots/final_plot_common_hwo.png')
    
        la.final_plot([hwo['st_spectype','sy_dist'],StarCat4['class_temp','dist_st_value'],
                  hpic['st_spectype','sy_dist']],
                  ['SPORES','LIFE-StarCat4','HPIC'],path='plots/final_plot_hwo_life.png')
    return
