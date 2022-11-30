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

	<table id="object" onDisk="True" adql="True">
                <meta name="title">Object table</meta>
                <meta name="description">
                A list of the object parameters.</meta>
		<primary>object_id</primary>
                <column name="object_id" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="object_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="type" type="text"
                        ucd="meta.record"
                        tablehead="type"
                        description="Object type."
                        required="True"
                        verbLevel="1"/>
        </table> 
       
        
        <data id="import_object_table">
        	<sources>data/objects.xml</sources>
		<!--	Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>
                <make table="object"> 
                	 <rowmaker idmaps="*"/>  
                </make>                           		
	</data>	 
	               
	<table id="star_basic" onDisk="True" adql="True" mixin="//scs#q3cindex">
		<meta name="title">Basic stellar parameters</meta>
		<meta name="description">
		A list of all basic stellar parameters.</meta>
		<primary>object_idref</primary>
		<foreignKey source="object_idref" inTable="object"
                        dest="object_id" /> 
		<!-- <foreignKey source="coord_source_idref" inTable="source"
                        dest="source_id" /> -->
                <!-- <foreignKey source="plx_source_idref" inTable="source"
                        dest="source_id" /> -->
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
                        verbLevel="1"/>
                <column name="plx_value" type="text"
                        ucd="pos.parallax"
                        tablehead="plx_value"
                        description="Parallax value."
                        verbLevel="1"/>
                <column name="plx_err" type="text"
                        ucd="pos.parallax"
                        tablehead="plx_err"
                        description="Parallax uncertainty."
                        verbLevel="1"/>
                <column name="plx_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="plx_source_id"
                        description="Source identifier corresponding 
                        to the parallax parameters."
                        verbLevel="1"/>
	</table>
	
	<data id="import_star_basic_table">
		<sources>data/stars.xml</sources>
		<!--Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>

                <make table="star_basic">
                        <rowmaker idmaps="*">
                        	<map dest="coord_ra" src="ra"/>
                        	<map dest="coord_source_idref"
                        	 src="coord_source_idref"
                        	 nullExpr="0" />
                        	<map dest="plx_source_idref"
                        	 src="plx_source_idref"
                        	 nullExpr="0" />
                        	<map dest="coord_dec" src="dec"/>
                        </rowmaker>
                </make>                                       		
	</data>	
	
	<table id="planet_basic" onDisk="True" adql="True">
		<meta name="title">Basic planetary parameters</meta>
		<meta name="description">
		A list of all basic planetary parameters.</meta>
		<primary>object_idref</primary>
		<foreignKey source="object_idref" inTable="object"
                        dest="object_id" /> 
		<!-- <foreignKey source="bestmass_source_idref" inTable="source"
                        dest="source_id" /> -->
                <column name="object_idref" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="object_idref"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
		<column name="bestmass" type="double precision"
			ucd="phys.mass" unit="Mjup" 
			tablehead="bestmass" 
			description="ExoMercat bestmass parameter." 
			verbLevel="1"/>
		<column name="bestmass_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="bestmass_source_idref"
                        description="Identifier of the source of the bestmass
                         parameter."
                        verbLevel="1"/>	
	</table>
    
	<data id="import_plant_basic_table">
		<sources>data/planets.xml</sources>
		<!--	Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>

		<make table="planet_basic">
                        <rowmaker idmaps="*">
                        	<map dest="bestmass_source_idref"
                        	 src="bestmass_source_idref"
                        	 nullExpr="0" />
                        </rowmaker>
                </make>       
                                       		
	</data>	

	<table id="disk_basic" onDisk="True" adql="True">
		<meta name="title">Basic disk parameters</meta>
		<meta name="description">
		A list of all basic disk parameters.</meta>
		<primary>object_idref</primary>
		<foreignKey source="object_idref" inTable="object"
                        dest="object_id" /> 
		<!-- <foreignKey source="ref" inTable="source"
                        dest="source_id" />  -->
                <column name="object_idref" type="text"
                        ucd="meta.id;meta.main"
                        tablehead="object_idref"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
		<column name="radius" type="double precision"
			ucd="phys.siye.radius"  
			tablehead="radius" 
			description="Black body radius." 
			verbLevel="1"/>
		<column name="radius_error" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="radius_error"
                        description="Uncertainty of the radius parameter."
                        verbLevel="1"/>
                <column name="_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="mass_source_idref"
                        description="Identifier of the source of the disk
                         parameters."
                        required="True"
                        verbLevel="1"/>	
	</table>
    
	<data id="import_disk_basic_table">
		<sources>data/disks.xml</sources>
		<!--	Data acquired through personal communication with 
		Grant Kennedy -->
		<voTableGrammar/>

		<make table="disk_basic">
                        <rowmaker idmaps="*">
                        	<map dest="radius">
                        		parseWithNull(@rdisk_bb,float,"None")</map>
                        	<map dest="radius_error"
                        	 src="e_rdisk_bb"/>
                        </rowmaker>
                </make>                                             		
	</data>	

	
	<table id="h_link" onDisk="True" adql="True">
                <meta name="title">Object relation table</meta>
                <meta name="description">
                A list of the object relations.</meta>
		<primary>parent_object_idref,child_object_idref,membership,h_link_source_idref</primary>
                <column name="parent_object_idref" type="text"
                        ucd="meta.id.parent;meta.main"
                        tablehead="parent"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="child_object_idref" type="text"
                        ucd="meta.id"
                        tablehead="child"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="membership" type="text"
                        ucd="meta.record"
                        tablehead="membership"
                        description="Membership probability."
                        verbLevel="1"/>
                <column name="h_link_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="h_link_source_idref"
                        description="Identifier of the source of the relationship
                         parameters."
                        required="True"
                        verbLevel="1"/>	
        </table> 
       
        
        <data id="import_h_link_table">
        	<sources>data/h_link.xml</sources>
		<!--	Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>
                <make table="h_link">  
                	 <rowmaker idmaps="*">
                	 </rowmaker>
                </make>                           		
	</data>
		 
	<table id="ident" onDisk="True" adql="True">
                <meta name="title">Object identifiers table</meta>
                <meta name="description">
                A list of the object identifiers.</meta>
		<primary>object_idref,id,id_source_idref</primary>
                <column name="object_idref" type="text"
                        ucd="meta.id"
                        tablehead="object_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="id" type="text"
                        ucd="meta.id"
                        tablehead="id"
                        description="Object identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="id_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="id_source_idref"
                        description="Identifier of the source of the identifier."
                        required="True"
                        verbLevel="1"/>	
        </table> 
       
        
        <data id="import_ident_table">
        	<sources>data/ident.xml</sources>
		<!--	Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>
                <make table="ident">  
                	 <rowmaker idmaps="*">
                	 </rowmaker>
                </make>                           		
	</data>		

	<table id="mesDist" onDisk="True" adql="True">
                <meta name="title">Distance measurement table</meta>
                <meta name="description">
                A list of the stellar distance measurements.</meta>
		<primary>object_idref,dist,mesDist_source_idref</primary>
                <column name="object_idref" type="text"
                        ucd="meta.id"
                        tablehead="object_id"
                        description="Object internal identifier."
                        required="True"
                        verbLevel="1"/>
                <column name="dist" type="text"
                        ucd="pos.distance"
                        tablehead="dist"
                        description="Object distance."
                        verbLevel="1"/>
                <column name="mesDist_source_idref" type="text"
                        ucd="meta.record;meta.id"
                        tablehead="mesDist_source_idref"
                        description="Identifier of the source of the distance measurement."
                        verbLevel="1"/>	
        </table> 
       
        
        <data id="import_mesDist_table">
        	<sources>data/mesDist.xml</sources>
		<!--	Data acquired using the data_acquisition.py skript -->
		<voTableGrammar/>
                <make table="mesDist">  
                	 <rowmaker idmaps="*">
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


