<resource schema="life">
    <meta name="creationDate">2021-12-03T09:33:44Z</meta>

<!-- Adjust the following text when you'd like to tone the beta
    warning text up or down.  And don't indent it to avoid spurious
    whitespace. -->
<macDef name="betawarning">Note that LIFE's target database is living
data.  The content – and to some extent even structure – of these
tables may change at any time without prior warning.
</macDef>

    <meta name="title">The LIFE Target Star Database LIFETD</meta>
    <meta name="description">
    The LIFE Target Star Database contains information useful
    for the planned LIFE mission (mid-ir, nulling
    interferometer in space). It characterizes possible
    target systems including information about stellar,
    planetary and disk properties. The data itself is mainly
    a collection from different other catalogs.

    \betawarning</meta>
    <meta name="subject">stars</meta>
    <meta name="subject">exoplanets</meta>
    <meta name="subject">circumstellar-disks</meta>
    <meta name="subject">astrometry</meta>

    <meta name="creator">Menti, F.; Quanz, S.; LIFE Collaboration</meta>
    <meta name="instrument">LIFE</meta> 
    <meta name="contentLevel">Research</meta>
    <meta name="type">Archive</meta>  

    <table id="source" onDisk="True" adql="True">
        <meta name="title">Source Table</meta>
        <meta name="description">
        A list of all the sources for the parameters in the other
        tables.

        \betawarning</meta>
        <primary>source_id</primary>
        <column name="source_id" type="integer"
            ucd="meta.id"
            tablehead="source_id"
            description="Source identifier."
            required="True"
            verbLevel="1"/>
        <column name="ref" type="text"
            ucd="meta.ref"
            tablehead="ref"
            description="Reference, bibcode if possible."
            verbLevel="1"/> 
        <column name="provider_name" type="text"
            ucd="meta.bib.author"
            tablehead="provider_name"
            description="Provider / parameter acquisition source."
            required="True"
            verbLevel="1"/>
        <column name="provider_url" type="text"
            ucd="meta.ref.url"
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
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="source">  
            <rowmaker idmaps="*">
            </rowmaker>
        </make>  
    </data>

    <table id="object" onDisk="True" adql="True">
        <meta name="title">Object table</meta>
        <meta name="description">
        A list of the astrophysical objects.

        \betawarning</meta>
        <primary>object_id</primary>
        <column name="object_id" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_id"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="type" type="text"
            ucd="src.class"
            tablehead="type"
            description="Object type (sy=multi-object system, st=star, 
            pl=planet and di=disk)."
            required="True"
            verbLevel="1"/>
        <column name="main_id" type="text"
            ucd="meta.id"
            tablehead="main_id"
            description="Main object identifier."
            required="True"
            verbLevel="1"/>
        <column name="ids" type="text"
            ucd="meta.id"
            tablehead="ids"
            description="All identifiers of the object separated by '|'."
            required="True"
            verbLevel="1"/>
    </table> 
    
    <data id="import_object_table">
        <sources>data/objects.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="object"> 
             <rowmaker idmaps="*"/>  
        </make>                       
    </data>     
               
    <table id="star_basic" onDisk="True" adql="True" mixin="//scs#q3cindex">
        <meta name="title">Basic stellar parameters</meta>
        <meta name="description">
        A list of all basic stellar parameters.

        \betawarning</meta>
        <primary>object_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />
  
        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="coo_ra" type="double precision"
            ucd="pos.eq.ra;meta.main" unit="deg" 
            tablehead="RA (ICRS)" 
            description="Right Ascension" 
            verbLevel="1"/>
        <column name="coo_dec" type="double precision" 
            ucd="pos.eq.dec;meta.main" unit="deg"
            tablehead="Dec (ICRS)" 
            description="Declination"
            verbLevel="1"/>
        <column name="coo_err_angle" type="smallint" 
            ucd="pos.posAng;pos.errorEllipse;pos.eq" unit="deg"
            tablehead="coo_err_angle" 
            description="Coordinate error angle"
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="coo_err_maj" type="real" 
            ucd="phys.angSize.smajAxis;pos.errorEllipse;pos.eq" 
            unit="mas"
            tablehead="coo_err_maj" 
            description="Coordinate error major axis"
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="coo_err_min" type="real" 
            ucd="phys.angSize.sminAxis;pos.errorEllipse;pos.eq" 
            unit="mas"
            tablehead="coo_err_min" 
            description="Coordinate error minor axis"
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="coo_qual" type="text" 
            ucd="meta.code.qual;pos.eq"
            tablehead="coo_qual" 
            description="Coordinate quality"
            verbLevel="1"/>
        <column name="coo_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="coord_source_id"
            description="Source identifier corresponding 
            to the position (coo) parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="plx_value" type="double precision"
            ucd="pos.parallax" unit="mas"
            tablehead="plx_value"
            description="Parallax value."
            verbLevel="1"/>
        <column name="plx_err" type="double precision"
            ucd="stat.error;pos.parallax" unit="mas"
            tablehead="plx_err"
            description="Parallax uncertainty."
            verbLevel="1"/>
        <column name="plx_qual" type="text" 
            ucd="meta.code.qual;pos.parallax"
            tablehead="plx_qual" 
            description="Parallax quality"
            verbLevel="1"/>
        <column name="plx_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="plx_source_id"
            description="Source identifier corresponding 
            to the parallax parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="dist_value" type="double precision"
            ucd="pos.distance"
            tablehead="dist"
            description="Object distance."
            required="True"
            verbLevel="1"/>
        <column name="dist_err" type="double precision"
            ucd="stat.error;pos.distance"
            tablehead="dist_err"
            description="Object distance error."
            verbLevel="1"/>
        <column name="dist_qual" type="text" 
            ucd="meta.code.qual;pos.distance"
            tablehead="dist_qual" 
            description="Distance quality"
            verbLevel="1"/>
        <column name="dist_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="dist_source_idref"
            description="Identifier of the source of the 
                distance parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>
    
    <data id="import_star_basic_table">
        <sources>data/star_basic.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="star_basic">
            <rowmaker idmaps="*">
                <map dest="coo_source_idref"
                 src="coo_source_idref"
                 nullExpr="0" />
                <map dest="plx_source_idref"
                 src="plx_source_idref"
                 nullExpr="0" />
            </rowmaker>
        </make>                               
    </data>    
    
    <table id="planet_basic" onDisk="True" adql="True">
        <meta name="title">Basic planetary parameters</meta>
        <meta name="description">
        A list of all basic planetary parameters.
        
        \betawarning
        </meta>
        <primary>object_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" /> 

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="mass_val" type="double precision"
            ucd="phys.mass" unit="'jupiterMass'" 
            tablehead="mass_val" 
            description="Mass" 
            verbLevel="1"/>
        <column name="mass_err" type="double precision"
            ucd="stat.error;phys.mass" unit="'jupiterMass'"
            tablehead="mass_err"
            description="Mass error"
            verbLevel="1"/>
        <column name="mass_qual" type="text" 
            ucd="meta.code.qual;phys.mass"
            tablehead="mass_qual" 
            description="Mass quality"
            verbLevel="1"/>
        <column name="mass_rel" type="text" 
            ucd="phys.mass;arith.ratio"
            tablehead="mass_rel" 
            description="Mass relation defining upper / lower limit or exact 
            measurement through '&lt;', '>', and '='."
            verbLevel="1"/>
        <column name="mass_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mass_source_idref"
            description="Identifier of the source of the mass
             parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>
    
    <data id="import_plant_basic_table">
        <sources>data/planet_basic.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="planet_basic">
            <rowmaker idmaps="*">
                <map dest="mass_source_idref"
                 src="mass_source_idref"
                 nullExpr="0" />
            </rowmaker>
        </make>                           
    </data>    

    <table id="disk_basic" onDisk="True" adql="True">
        <meta name="title">Basic disk parameters</meta>
        <meta name="description">
        A list of all basic disk parameters.
        
        \betawarning
        </meta>
        <primary>object_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" /> 
        
        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="rad_value" type="double precision"
            ucd="phys.size.radius"  unit="AU"
            tablehead="radius" 
            description="Black body radius." 
            verbLevel="1"/>
        <column name="rad_err" type="double precision"
            ucd="stat.error;phys.size.radius" unit="AU"
            tablehead="radius_error"
            description="Radius error"
            verbLevel="1"/>
        <column name="rad_qual" type="text" 
            ucd="meta.code.qual;phys.size.radius"
            tablehead="rad_qual" 
            description="Radius quality"
            verbLevel="1"/>
        <column name="rad_rel" type="text" 
            ucd="phys.size.radius;arith.ratio"
            tablehead="rad_rel" 
            description="Radius relation defining upper / lower limit or exact 
            measurement through '&lt;' ,'>', and '='."
            verbLevel="1"/>
        <column name="rad_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mass_source_idref"
            description="Identifier of the source of the disk
                parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>
    
    <data id="import_disk_basic_table">
        <sources>data/disk_basic.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="disk_basic">
            <rowmaker idmaps="*">
                <map dest="rad_value">
                    parseWithNull(@rad_value,float,"None")
                </map>
            </rowmaker>
        </make>                                 
    </data>    

    
    <table id="h_link" onDisk="True" adql="True">
        <meta name="title">Object relation table</meta>
        <meta name="description">
        This table links subordinate objects (e.g. a planets of a star, or
        a star in a multiple star system) to their parent objects.

        \betawarning
        </meta>
        <primary>parent_object_idref,child_object_idref,
            h_link_source_idref
        </primary>
        <column name="parent_object_idref" type="integer"
            ucd="meta.id.parent;meta.main"
            tablehead="parent"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="child_object_idref" type="integer"
            ucd="meta.id"
            tablehead="child"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="membership" type="smallint"
            ucd="meta.record"
            tablehead="membership"
            description="Membership probability."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>      
        <column name="h_link_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="h_link_source_idref"
            description="Identifier of the source of the 
                relationship parameters."
            required="True"
            verbLevel="1"/>    
    </table> 
       
    <data id="import_h_link_table">
        <sources>data/h_link.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="h_link">  
             <rowmaker idmaps="*">
                 <map dest="membership"
                 src="membership"
                 nullExpr="-1" />
             </rowmaker>
        </make>                       
    </data>
         
    <table id="ident" onDisk="True" adql="True">
        <meta name="title">Object identifiers table</meta>
        <meta name="description">
        A list of the object identifiers.
        
        \betawarning
        </meta>
        <primary>object_idref,id,id_source_idref</primary>
        <column name="object_idref" type="integer"
            ucd="meta.id"
            tablehead="object_id"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="id" type="text"
            ucd="meta.id"
            tablehead="id"
            description="Object identifier."
            required="True"
            verbLevel="1"/>
        <column name="id_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="id_source_idref"
            description="Identifier of the source of the 
                identifier parameter."
            required="True"
            verbLevel="1"/>    
    </table> 
       
    <data id="import_ident_table">
        <sources>data/ident.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="ident">  
             <rowmaker idmaps="*">
             </rowmaker>
        </make>                       
    </data>        

    <table id="mes_dist" onDisk="True" adql="True">
        <meta name="title">Distance measurement table</meta>
        <meta name="description">
        A list of the stellar distance measurements.

        \betawarning
        </meta>
        <primary>object_idref,dist_value,dist_err,
            dist_source_idref
        </primary>
        
        <column name="object_idref" type="integer"
            ucd="meta.id"
            tablehead="object_id"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="dist_value" type="double precision"
            ucd="pos.distance" unit="AU"
            tablehead="dist_value"
            description="Object distance."
            required="True"
            verbLevel="1"/>
        <column name="dist_err" type="double precision"
            ucd="stat.error;pos.distance" unit="AU"
            tablehead="dist_err"
            description="Distance error."
            verbLevel="1"/>
        <column name="dist_qual" type="text" 
            ucd="meta.code.qual;pos.distance"
            tablehead="dist_qual" 
            description="Distance quality"
            verbLevel="1"/>
        <column name="dist_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mes_dist_source_idref"
            description="Identifier of the source of the 
                distance parameter."
            required="True"
            verbLevel="1"/>    
    </table> 
       
    <data id="import_mes_dist_table">
        <sources>data/mesDist.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="mes_dist">  
             <rowmaker idmaps="*">
             </rowmaker>
        </make>                       
    </data>                      
  
    <table id="mes_mass" onDisk="True" adql="True">
        <meta name="title">Mass measurement table</meta>
        <meta name="description">
        A list of the planetary mass measurements.

        \betawarning
        </meta>
        <primary>object_idref,mass_val,mass_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" /> 
            
        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the 
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="mass_val" type="double precision"
            ucd="phys.mass" unit="'jupiterMass'" 
            tablehead="mass_val" 
            description="Mass" 
            verbLevel="1"/>
        <column name="mass_err" type="double precision"
            ucd="stat.error;phys.mass" unit="'jupiterMass'"
            tablehead="mass_err"
            description="Mass error"
            verbLevel="1"/>
        <column name="mass_qual" type="text" 
            ucd="meta.code.qual;phys.mass"
            tablehead="mass_qual" 
            description="Mass quality"
            verbLevel="1"/>
        <column name="mass_rel" type="text" 
            ucd="phys.mass;arith.ratio"
            tablehead="mass_rel" 
            description="Mass relation defining upper / lower limit or exact 
            measurement through '&lt;', '>', and '='."
            verbLevel="1"/>
        <column name="mass_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mass_source_idref"
            description="Identifier of the source of the mass
             parameter."
            required="True"
            verbLevel="1"/>
    </table> 
       
    <data id="import_mes_mass_table">
        <sources>data/mesMass.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_mass">  
             <rowmaker idmaps="*">
             </rowmaker>
        </make>                       
    </data>                      
 
    <service id="cone" allowed="form,scs.xml">
        <meta name="shortName">lifetd cone</meta>
        <meta>
            testQuery.ra: 312.27
            testQuery.dec: 37.47
            testQuery.sr: 0.01
        </meta>
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
<!-- vim:et:sw=4:sta
-->
