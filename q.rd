<resource schema="lifetest">
	<meta name="creationDate">2021-12-03T09:33:44Z</meta>

	<meta name="title">The LIFE Target Star Database LIFETD</meta>
	<meta name="description">
	The LIFE Target Star Database contains information useful
	for the planned LIFE mission (mid-ir, nulling
	interferometer in space). It characterizes possible
	target systems including information about stellar,
	planetary and disk properties. The data itself is mainly
	a collection from different other catalogs.</meta>
	<meta name="subject">stars</meta>
	<meta name="subject">planets</meta>
	<meta name="subject">disks</meta>
	<meta name="subject">astrometry</meta>
	<meta name="subject">proper-motions</meta>

	<meta name="creator">Menti, F.; Quanz, S.; LIFE team</meta>
	<meta name="instrument">LIFE</meta> 
	<meta name="facility">ETH Zurich</meta>
	<meta name="source">LIFE source</meta>
	<meta name="contentLevel">Research</meta>
	<meta name="type">Archive</meta>  
	<meta name="coverage.waveband">Infrared</meta>  
  
	<table id="star_basic" onDisk="True" adql="True">
		<meta name="title">Basic stellar parameters</meta>
		<meta name="description">
		A list of all basic stellar parameters.</meta>
		<primary>star_id</primary>
                <column name="star_id" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="star_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
		<column name="ra" type="double precision"
			ucd="pos.eq.ra;meta.main" unit="deg" 
			tablehead="RA (ICRS)" 
			description="Right Ascension" 
			verbLevel="1">
		</column>
		<column name="dec" type="double precision" 
			ucd="pos.eq.dec;meta.main" unit="deg"
			tablehead="Dec (ICRS)" 
			description="Declination"
			verbLevel="1">
		</column>
		<column name="new_id" type="text"
			ucd="meta.record;meta.id"
                        tablehead="new_id"
                        description="New internal identifier."
                        required="True"
                	verbLevel="1"/>
			
	</table>
	<!-- <table id="planet_basic" onDisk="True" adql="True">
		<meta name="title">Basic planetarz parameters</meta>
		<meta name="description">
		A list of all basic planetarz parameters.</meta>
		<primary>planet_id</primary>
                <column name="planet_id" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="star_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
		<column name="mass" type="double precision"
			ucd="phys.mass" unit="Mjup" 
			tablehead="mass" 
			description="ExoMercat Bestmass parameter." 
			verbLevel="1">
		</column>
	</table> -->
	<table id="object" onDisk="True" adql="True">
                <meta name="title">Object table</meta>
                <meta name="description">
                A list of the object parameters.</meta>
                <primary>main_id</primary>
                <foreignKey dest="star_id" inTable="star_basic"
                        source="object_id" /> 
		<column name="object_id" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="object_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="main_id" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="main_id"
                        description="Main identifier."
                        required="True"
                        verbLevel="1"/>
        </table>
    
	<data id="import">
		<sources>
			<pattern>data/simbad_main_id.fits</pattern>
				<!--Data queried from SIMBAD:
					SELECT TOP 10 main_id,ra,dec,oid
					FROM basic
					WHERE basic.plx_value >=50. -->
			<!-- <pattern>data/exomercat_10_id.fits</pattern>
				Data queried from exomercat:
					SELECT TOP 10 id, name, bestmass
					FROM exomercat.exomercat  -->
		</sources>
		<fitsTableGrammar/>
		<!-- csvGrammar/ -->

                <make table="star_basic">
                        <rowmaker idmaps="*">
                        	<var key="new_id">"STAR%"%(@oid)</var>
                        	<map dest="star_id" src="oid"/>
                        </rowmaker>
                </make>
                <!--<make table="planet_basic">
                        <rowmaker idmaps="*">
                        	<map dest="planet_id" src="id"/>
                        	<map dest="mass" src="bestmass"/>
                        </rowmaker>
                </make> -->            
                <make table="object">
                        <rowmaker idmaps="*">
                        	<map dest="object_id" src="oid"/>
                        	<!-- <map dest="object_id" src="id"/>
                        	<map dest="main_id" src="name"/>  -->              
                        </rowmaker>
                </make>                       		
	</data>	
                                      
  
	<service id="cone" allowed="form,scs.xml">
		<meta name="shortName">service short name</meta>
		<publish render="form" sets="ivo_managed, local"/>
		<publish render="scs.xml" sets="ivo_managed"/>
		<scsCore queriedTable="star_basic">
		</scsCore>
	</service> 
  
	<regSuite title="life regression">
		<regTest title="LIFE form service appears to work."
			url="scs/form">
			<code>
				assert False
			</code>
		</regTest>
	</regSuite>
</resource>
