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
                <column name="ref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="ref"
                        description="Reference, bibcode if possible."
                        verbLevel="1"/> 
                <column name="provider_name" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="provider_name"
                        description="By what source the parameter was acquired."
                        required="True"
                        verbLevel="1"/>
                <column name="provider_url" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="provider_url"
                        description="Provider URL."
                        required="True"
                        verbLevel="1"/>
                <column name="provider_bibcode" type="text"
                        ucd="meta.bib.bibcode"
                        tablehead="provider_bibcode"
                        description="Provider bibcode."
                        required="True"
                        verbLevel="1"/>              
        </table>
        
        <data id="import_source_table">
        	<sources>data/sources.xml</sources>
		<!-- Data queried with data_acquisition.py -->
		<voTableGrammar/>
                <make table="source">  
                	<rowmaker idmaps="*">
                	</rowmaker>
                </make>  
        </data>
               
	<table id="star_basic" onDisk="True" adql="True" mixin="//scs#q3cindex">
		<meta name="title">Basic stellar parameters</meta>
		<meta name="description">
		A list of all basic stellar parameters.</meta>
		<primary>object_idref</primary>
		<!-- <foreignKey source="object_idref" inTable="object"
                        dest="object_id" /> -->
		<foreignKey source="coord_source_idref" inTable="source"
                        dest="source_id" />
                <foreignKey source="plx_source_idref" inTable="source"
                        dest="source_id" />
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
                        tablehead="coord_source_id"
                        description="Source identifier corresponding 
                        to the position (coord) parameters."
                        required="True"
                        verbLevel="1"/>
                <column name="plx_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="plx_source_id"
                        description="Source identifier corresponding 
                        to the parallax parameters."
                        required="True"
                        verbLevel="1"/>
	</table>
	
	<data id="import_star_basic_table">
		<sources>data/stars.xml</sources>
		<!--Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>

                <make table="star_basic">
                        <rowmaker idmaps="*">
                        	<var key="object_idref">"S%s"%@oid</var>
                        	<map dest="coord_ra" src="ra"/>
                        	<map dest="coord_dec" src="dec"/>
                        </rowmaker>
                </make>                                       		
	</data>	
	
	<table id="planet_basic" onDisk="True" adql="True">
		<meta name="title">Basic planetary parameters</meta>
		<meta name="description">
		A list of all basic planetary parameters.</meta>
		<primary>object_idref</primary>
		<foreignKey source="mass_source_idref" inTable="source"
                        dest="source_id" />
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
                        tablehead="mass_source_idref"
                        description="Identifier of the source of the mass
                         parameter."
                        required="True"
                        verbLevel="1"/>	
	</table>
    
	<data id="import_plant_basic_table">
		<sources>data/planets.xml</sources>
		<!--	Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>

		<make table="planet_basic">
                        <rowmaker idmaps="*">
                        	<map dest="mass" src="bestmass"/>
                        	<map dest="mass_source_idref"
                        	 src="bestmass_source_idref"/>
                        	<var key="object_idref">"P%s"%@id</var> 
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
