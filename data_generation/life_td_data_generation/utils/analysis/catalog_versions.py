import numpy as np #arrays
import astropy as ap #votables
import importlib #reloading external functions after modification
import seaborn as sns
import sys
from matplotlib import pyplot as plt

sys.path.append('../../implementation/life/data_generation/life_td_data_generation')

#self created modules
from utils.io import save, load, stringtoobject, Path
from utils.analysis.catalog_comparison import create_common
import utils.analysis.analysis as la
importlib.reload(la)#reload module after changing it
import starcat4 as ltc4
#import LIFE_StarCat4 as ltc4
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

def threecatboxplot(data,para,labels):
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
    bp = ax.boxplot(data, labels=labels)
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

def compare():
    para4=['coo_ra', 'coo_dec', 'plx_value', 'dist_st_value', 'coo_gal_l', 'coo_gal_b',
       'teff_st_value', 'radius_st_value','mass_st_value','sep_phys_value',  'mag_i_value', 'mag_j_value']
    paras=[para4,para4]
    labels=['StarCat5','StarCat4']
    paths=['catalogs/StarCat5','StarCat4']
    ltc_compare(paths, labels, paras)

def ltc_compare(paths, labels, paras):
    """
    This function compares StarCat4 to StarCat3 and HWO catalog SPORES.

    :param bool load: Defaults to True, determines if tables are created or loaded.
    :param float distance_cut_in_pc: Defaults to 30., determines distance cut of StarCat4 if load=False.
    """

    boxplots = True
    finalplots = True

    catalogs=[]
    for path in paths:
        catalogs.append(load([path],location=Path().additional_data)[0])

    common_0_1, common_1_0 = create_common(catalogs[0], catalogs[1])
    join_0_1 = ap.table.join(common_0_1, common_1_0, keys='main_id')
    save([join_0_1], [f'join_{labels[0]}_{labels[1]}'],
         location=Path().additional_data)


    if boxplots:
        #comparing only objects that are in both catalogs
        df_old_new=prepare_dataframe_for_sns(join_0_1,[paras[0],paras[1]],
                                             [labels[0],labels[1]])

        for y in paras[0]:
            snsplot(df_old_new,y,f'plots/sns_{labels[0]}_{labels[1]}{y}')

        y_axis=['Stellar_Distance','Stellar_Mass',
                'Stellar_Radius','Stellar_Effective_Temperature']
        for i in range(len(paras[0])):
        #comparing all objects in the catalogs -> wait, I want only two cats -> add very old
            data = [catalogs[0][paras[0][i]],catalogs[1][paras[1][i]]]
            threecatboxplot(data,y_axis[i],labels)

    if finalplots:

        la.final_plot([new['class_temp','dist_st_value'],old['sim_sptype','distance']],
                  ['LIFE-StarCat4','LIFE-StarCat3'],path='plots/final_plot_34.png')

        la.final_plot([common_new_and_old['class_temp','dist_st_value'],common_old_and_new['sim_sptype','distance']],
                  ['common_LIFE-StarCat4','common_LIFE-StarCat3'],path='plots/final_plot_common_34.png')

    return
