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
       
	<table id="star_basic" onDisk="True" adql="True" mixin="//scs#q3cindex">
		<meta name="title">Basic stellar parameters</meta>
		<meta name="description">
		A list of all basic stellar parameters.</meta>
		<primary>star_id</primary>
		<!-- <foreignKey source="object_idref" inTable="object"
                        dest="object_id" /> -->
                <column name="star_id" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="star_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="object_idref" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="object_idref"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
		<column name="coord_ra" type="double precision"
			ucd="pos.eq.ra;meta.main" unit="deg" 
			tablehead="RA (ICRS)" 
			description="Right Ascension" 
			verbLevel="1"/>
		<column name="coord_dec" type="double precision" 
			ucd="pos.eq.dec;meta.main" unit="deg"
			tablehead="Dec (ICRS)" 
			description="Declination"
			verbLevel="1"/>
		<column name="coord_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="coord_source_idref"
                        description="Identifier of the source of the coordinate
                         parameter."
                        required="True"
                        verbLevel="1"/>	
	</table>
	
	<data id="import_star_basic_table">
		<sources>data/star.fits</sources>
		<!--Data queried from SIMBAD:
				SELECT TOP 10 main_id,ra,dec,oid,
				 coo_bibcode
				FROM basic
				WHERE basic.plx_value >=50. -->
		<fitsTableGrammar/>

                <make table="star_basic">
                        <rowmaker idmaps="*">
                        	<var key="object_idref">"STAR-%s"%@oid</var>
                        	<map dest="star_id" src="oid"/>
                        	<map dest="coord_ra" src="ra"/>
                        	<map dest="coord_dec" src="dec"/>
                        	<map dest="coord_source_idref" 
                        		src="coo_bibcode"/>
                        </rowmaker>
                </make>                                       		
	</data>	
	
	<table id="planet_basic" onDisk="True" adql="True">
		<meta name="title">Basic planetary parameters</meta>
		<meta name="description">
		A list of all basic planetary parameters.</meta>
		<primary>planet_id</primary>
                <column name="planet_id" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="planet_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="object_idref" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="object_idref"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
		<column name="mass" type="double precision"
			ucd="phys.mass" unit="Mjup" 
			tablehead="mass" 
			description="ExoMercat Bestmass parameter." 
			verbLevel="1"/>
		<column name="mass_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="coord_source_idref"
                        description="Identifier of the source of the mass
                         parameter."
                        required="True"
                        verbLevel="1"/>	
	</table>
    
	<data id="import_plant_basic_table">
		<sources>
			<pattern>data/planets.fits</pattern>
			<!--	Data queried from exomercat:
					SELECT TOP 20 id, name, bestmass, bestmass_url
					FROM exomercat.exomercat  -->
		</sources>
		<fitsTableGrammar/>

		<make table="planet_basic">
                        <rowmaker idmaps="*">
                        	<map dest="planet_id" src="id"/>
                        	<map dest="mass" src="bestmass"/>
                        	<map dest="mass_source_idref" src="bestmass_url"/>
                        	<var key="object_idref">"PLANET-%s"%@id</var> 
                        </rowmaker>
                </make>       
                                       		
	</data>	


	<table id="object" onDisk="True" adql="True">
                <meta name="title">Object table</meta>
                <meta name="description">
                A list of the object parameters.</meta>
		<primary>object_id</primary>
                <column original="star_basic.object_idref" name="object_id"/>        	
        	
        	<viewStatement>
        		CREATE VIEW \curtable AS ( 
        			SELECT \colNames FROM
        				(SELECT object_idref AS object_id
        					FROM \schema.star_basic) as subquery1
        				UNION
        				(SELECT object_idref as object_id
        					FROM \schema.planet_basic))
        	</viewStatement> 
        </table> 
       
        
        <data id="import_object_table">
                <make table="object">      
                </make>                           		
	</data>	 
	
	<table id="source" onDisk="True" adql="True">
		<meta name="title">Source Table</meta>
		<meta name="description">
		A list of all the sources for the parameters in the other
		 tables.</meta>
		<primary>source_id</primary>
                <column name="source_id" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="source_id"
                        description="Source identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="bibcode" type="text"
                        ucd="meta.bib.bibcode"
                        tablehead="bibcode"
                        description="Bibcode"
                        verbLevel="1"/> 
                <column name="acquired_through" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="acquired_through"
                        description="By what source the parameter was acquired."
                        required="True"
                        verbLevel="1"/>
                         
        </table>
        
        <data id="import_source_table">
        	<sources>data/source.fits</sources>
		<!-- Data queried from SIMBAD:
				SELECT DISTINCT * FROM (
					SELECT TOP 10 coo_bibcode
					FROM basic
					WHERE basic.plx_value >=50.) AS subquery
				JOIN ref ON bibcode=coo_bibcode -->
		<fitsTableGrammar/>

                <make table="source">  
                	<rowmaker idmaps="*">
                		<var key="acquired_through">"simbad"</var>  
                		<map dest="source_id" src="coo_bibcode"/>
                		<map dest="bibcode" src="coo_bibcode"/>
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
