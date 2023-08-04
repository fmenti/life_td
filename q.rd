<resource schema="life_td" resdir=".">
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

    <meta name="doi">10.21938/ke_e6lzO_jjX_vzvVIcwZA</meta>

    <meta name="creator">Menti, F.; Quanz, S.; Alei E.; Caballero J.;
    Demleitner M.; Garcia Munoz A.; Kennedy G.; Schmitt U.; Stassun K.; 
    Wyatt M.; LIFE Collaboration</meta>
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
            description="Name for the service through which the data was
            acquired. SIMBAD: service =
            http://simbad.u-strasbg.fr:80/simbad/sim-tap,
            bibcode = 2000A&amp;AS..143....9W. ExoMercat: service =
            http://archives.ia2.inaf.it/vo/tap/projects,
            bibcode = 2020A&amp;C....3100370A. Everything else is acquired
            through private communication."
            required="True"
            verbLevel="1"/>
    </table>

    <data id="import_source">
        <sources>data/sources.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="source">
            <rowmaker idmaps="*">
                <map key="ref" nullExpr="'nan'"/>
                <map key="source_id" nullExpr="0"/>
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

    <data id="import_object">
        <sources>data/objects.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="object">
             <rowmaker idmaps="*"/>
        </make>
    </data>
    
    <table id="provider" onDisk="True" adql="True">
        <meta name="title">Provider Table</meta>
        <meta name="description">
        A list of all the providers for the parameters in the other
        tables.

        \betawarning</meta>
        <primary>provider_name</primary>
        <column name="provider_name" type="text"
            ucd="meta.bib.author"
            tablehead="provider_name"
            description="Name for the service through which the data was
            acquired."
            required="True"
            verbLevel="1"/>
        <column name="provider_url" type="text"
            ucd="meta.ref.url"
            tablehead="provider_url"
            description="Service through which the data was
            acquired."
            verbLevel="1"/>
        <column name="provider_bibcode" type="text"
            ucd="meta.ref"
            tablehead="provider_bibcode"
            description="Reference, bibcode if possible."
            verbLevel="1"/>
        <column name="provider_access" type="text"
            ucd="time"
            tablehead="provider_access"
            description="Date of access to provider."
            verbLevel="1"/>
    </table>

    <data id="import_provider">
        <sources>data/provider.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="provider">
            <rowmaker idmaps="*">
                <map key="provider_url" nullExpr="'None'"/>
                <map key="provider_bibcode" nullExpr="'None'"/>
            </rowmaker>
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
            verbLevel="1" displayHint="sf=7"/>
        <column name="coo_dec" type="double precision"
            ucd="pos.eq.dec;meta.main" unit="deg"
            tablehead="Dec (ICRS)"
            description="Declination"
            verbLevel="1" displayHint="sf=7"/>
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
            description="Coordinate quality (A:best, E:worst)"
            verbLevel="1"/>
        <column name="coo_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="coo_source_id"
            description="Source identifier corresponding
            to the position (coo) parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="coo_gal_l" type="double precision"
            ucd="pos.galactic.lon" unit="deg"
            tablehead="coo_gal_l"
            description="Galactical longitude."/>
        <column name="coo_gal_b" type="double precision"
            ucd="pos.galactic.lat" unit="deg"
            tablehead="coo_gal_b"
            description="Galactical latitude."/>
        <column name="coo_gal_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="coo_gal_source_id"
            description="Source identifier corresponding
            to the position (coo_gal) parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="mag_i_value" type="double precision"
            ucd="phot.mag;em.opt.I" unit=""
            tablehead="mag_i_value"
            description="Magnitude in I filter."/>
        <column name="mag_i_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mag_i_source_id"
            description="Source identifier corresponding
            to the Magnitude in I filter parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="mag_j_value" type="double precision"
            ucd="phot.mag;em.IR.J" unit=""
            tablehead="mag_j_value"
            description="Magnitude in J filter."/>
        <column name="mag_j_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mag_j_source_id"
            description="Source identifier corresponding
            to the Magnitude in J filter parameters."
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
            description="Parallax quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="plx_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="plx_source_id"
            description="Source identifier corresponding
            to the parallax parameters."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="dist_st_value" type="double precision"
            ucd="pos.distance" unit="pc"
            tablehead="Dist"
            description="Object distance."
            verbLevel="1" displayHint="sf=2"/>
        <column name="dist_st_err" type="double precision"
            ucd="stat.error;pos.distance" unit="pc"
            tablehead="dist_st_err"
            description="Object distance error."
            verbLevel="1"/>
        <column name="dist_st_qual" type="text"
            ucd="meta.code.qual;pos.distance"
            tablehead="dist_st_qual"
            description="Distance quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="dist_st_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="dist_st_source_idref"
            description="Identifier of the source of the
                distance parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="sptype_string" type="text"
            ucd="src.spType" 
            tablehead="sptype"
            description="Object spectral type MK."
            verbLevel="1" displayHint="sf=2"/>
        <column name="sptype_err" type="double precision"
            ucd="stat.error;src.spType"
            tablehead="sptype_err"
            description="Object spectral type MK error."
            verbLevel="1"/>
        <column name="sptype_qual" type="text"
            ucd="meta.code.qual;src.spType"
            tablehead="sptype_qual"
            description="Spectral type MK quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="sptype_source_idref" type="integer"
            ucd="meta.ref;src.spType"
            tablehead="sptype_source_idref"
            description="Identifier of the source of the
                spectral type MK parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="class_temp" type="text"
            ucd="src.spType" 
            tablehead="temp class"
            description="Object spectral type MK temperature class."
            verbLevel="1" displayHint="sf=2"/>
        <column name="class_temp_nr" type="text"
            ucd="src.spType" 
            tablehead="temp class nr"
            description="Object spectral type MK temperature class number."
            verbLevel="1" displayHint="sf=2"/>
        <column name="class_lum" type="text"
            ucd="src.spType" 
            tablehead="lum class"
            description="Object spectral type MK luminosity class."
            verbLevel="1" displayHint="sf=2"/>
        <column name="class_source_idref" type="integer"
            ucd="meta.ref;src.spType"
            tablehead="class_source_idref"
            description="Identifier of the source of the
                spectral type MK classification parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="teff_st_value" type="double precision"
            ucd="phys.temperature.effective" unit="K"
            tablehead="teff_st"
            description="Object effective temperature."
            verbLevel="1" displayHint="sf=2"/>
        <column name="teff_st_err" type="double precision"
            ucd="stat.error;phys.temperature.effective" unit="K"
            tablehead="teff_st_err"
            description="Object effective temperature error."
            verbLevel="1"/>
        <column name="teff_st_qual" type="text"
            ucd="meta.code.qual;phys.temperature.effective"
            tablehead="teff_st_qual"
            description="Effective temperature quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="teff_st_source_idref" type="integer"
            ucd="meta.ref;phys.temperature.effective"
            tablehead="teff_st_source_idref"
            description="Identifier of the source of the
                effective temperature parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="radius_st_value" type="double precision"
            ucd="phys.size.radius" unit="solRad"
            tablehead="radius_st"
            description="Object radius."
            verbLevel="1" displayHint="sf=2"/>
        <column name="radius_st_err" type="double precision"
            ucd="stat.error;phys.size.radius" unit="solRad"
            tablehead="radius_st_err"
            description="Object radius error."
            verbLevel="1"/>
        <column name="radius_st_qual" type="text"
            ucd="meta.code.qual;phys.size.radius"
            tablehead="radius_st_qual"
            description="Radius quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="radius_st_source_idref" type="integer"
            ucd="meta.ref;phys.size.radius"
            tablehead="radius_st_source_idref"
            description="Identifier of the source of the
                radius parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="mass_st_value" type="double precision"
            ucd="phys.mass" unit="solMass"
            tablehead="mass_st"
            description="Object mass."
            verbLevel="1" displayHint="sf=2"/>
        <column name="mass_st_err" type="double precision"
            ucd="stat.error;phys.mass" unit="solMass"
            tablehead="mass_st_err"
            description="Object mass error."
            verbLevel="1"/>
        <column name="mass_st_qual" type="text"
            ucd="meta.code.qual;phys.mass"
            tablehead="mass_st_qual"
            description="Mass quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="mass_st_source_idref" type="integer"
            ucd="meta.ref;phys.mass"
            tablehead="mass_st_source_idref"
            description="Identifier of the source of the
                mass parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="binary_flag" type="text"
            ucd="meta.code.multip" 
            tablehead="binary_flag"
            description="Binary flag."
            verbLevel="1" displayHint="sf=2"/>
        <column name="binary_qual" type="text"
            ucd="meta.code.qual"
            tablehead="binary_qual"
            description="Binary quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="binary_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="binary_source_idref"
            description="Identifier of the source of the
                binary_flag parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
        <column name="sep_phys_value" type="double precision"
            ucd="pos.angDistance" unit="arcsec"
            tablehead="Phys. separation"
            description="Phys. separation of binary."
            verbLevel="1" displayHint="sf=2"/>
        <column name="sep_phys_err" type="double precision"
            ucd="stat.error;pos.angDistance" unit="arcsec"
            tablehead="sep_phys_err"
            description="Object phys. separation error."
            verbLevel="1"/>
        <column name="sep_phys_qual" type="text"
            ucd="meta.code.qual;pos.angDistance"
            tablehead="sep_phys_qual"
            description="Phys. separation quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="sep_phys_source_idref" type="integer"
            ucd="meta.ref;pos.angDistance"
            tablehead="sep_phys_source_idref"
            description="Identifier of the source of the
                sep_phys parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>

    <data id="import_star_basic">
        <sources>data/star_basic.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="star_basic">
            <rowmaker idmaps="*">
                <map key="coo_ra" nullExpr="1e+20" />
                <map key="coo_dec" nullExpr="1e+20" />
                <map key="coo_err_angle" nullExpr="1e+20" />
                <map key="coo_err_maj" nullExpr="1e+20" />
                <map key="coo_err_min" nullExpr="1e+20" />
                <map key="coo_qual" nullExpr="'?'" />
                <map key="coo_source_idref" nullExpr="999999" />
                <map key="coo_gal_l" nullExpr="1e+20" />
                <map key="coo_gal_b" nullExpr="1e+20" />
                <!-- <map key="coo_gal_err_angle" nullExpr="1e+20" />
                <map key="coo_gal_err_maj" nullExpr="1e+20" />
                <map key="coo_gal_err_min" nullExpr="1e+20" />
                <map key="coo_gal_qual" nullExpr="'?'" /> -->
                <map key="coo_gal_source_idref" nullExpr="999999" />
                <map key="plx_value" nullExpr="1e+20" />
                <map key="plx_err" nullExpr="1e+20" />
                <map key="plx_qual" nullExpr="'?'" />
                <map key="plx_source_idref" nullExpr="999999" />
                <map key="mag_i_value" nullExpr="1e+20" />
                <map key="mag_i_source_idref" nullExpr="999999" />
                <map key="mag_j_value" nullExpr="1e+20" />
                <map key="mag_j_source_idref" nullExpr="999999" />
                <map key="dist_st_value" nullExpr="1e+20" />
                <map key="dist_st_err" nullExpr="1e+20" />
                <map key="dist_st_qual" nullExpr="'?'" />
                <map key="dist_st_source_idref" nullExpr="999999" />
                <map key="sptype_string" nullExpr="'None'" />
                <map key="sptype_err" nullExpr="1e+20" />
                <map key="sptype_qual" nullExpr="'?'" />
                <map key="sptype_source_idref" nullExpr="999999" />
                <map key="class_temp" nullExpr="'?'" />
                <map key="class_temp_nr" nullExpr="'?'" />
                <map key="class_lum" nullExpr="'?'" />
                <map key="class_source_idref" nullExpr="999999" />
                <map key="teff_st_value" nullExpr="1e+20" />
                <map key="teff_st_err" nullExpr="1e+20" />
                <map key="teff_st_qual" nullExpr="'?'" />
                <map key="teff_st_source_idref" nullExpr="999999" />
                <map key="radius_st_value" nullExpr="1e+20" />
                <map key="radius_st_err" nullExpr="1e+20" />
                <map key="radius_st_qual" nullExpr="'?'" />
                <map key="radius_st_source_idref" nullExpr="999999" />
                <map key="mass_st_value" nullExpr="1e+20" />
                <map key="mass_st_err" nullExpr="1e+20" />
                <map key="mass_st_qual" nullExpr="'?'" />
                <map key="mass_st_source_idref" nullExpr="999999" />
                <map key="binary_flag" nullExpr="'None'" />
                <map key="binary_qual" nullExpr="'?'" />
                <map key="binary_source_idref" nullExpr="999999" />
                <map key="sep_phys_value" nullExpr="1e+20" />
                <map key="sep_phys_err" nullExpr="1e+20" />
                <map key="sep_phys_qual" nullExpr="'?'" />
                <map key="sep_phys_source_idref" nullExpr="999999" />
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
        <column name="mass_pl_value" type="double precision"
            ucd="phys.mass" unit="'jupiterMass'"
            tablehead="Planet Mass"
            description="Mass"
            verbLevel="1" displayHint="sf=3"/>
        <column name="mass_pl_err" type="double precision"
            ucd="stat.error;phys.mass" unit="'jupiterMass'"
            tablehead="Err. Mass"
            description="Mass error"
            verbLevel="1" displayHint="sf=3"/>
        <column name="mass_pl_qual" type="text"
            ucd="meta.code.qual;phys.mass"
            tablehead="Quality Mass"
            description="Mass quality (A:best, E:worst)"
            verbLevel="1"/>
        <column name="mass_pl_rel" type="text"
            ucd="phys.mass;arith.ratio"
            tablehead="mass_rel"
            description="Mass relation defining upper / lower limit or exact
            measurement through '&lt;', '>', and '='."
            verbLevel="1"/>
        <column name="mass_pl_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mass_pl_source_idref"
            description="Identifier of the source of the mass
             parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>

    <data id="import_planet_basic">
        <sources>data/planet_basic.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="planet_basic">
            <rowmaker idmaps="*">
                <map key="mass_pl_value" nullExpr="1e+20" />
                <map key="mass_pl_err" nullExpr="1e+20" />
                <map key="mass_pl_qual" nullExpr="'?'" />
                <map key="mass_pl_source_idref" nullExpr="999999" />
                <map key="mass_pl_rel" nullExpr="'?'"/>
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
            description="Radius quality (A:best, E:worst)"
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

    <data id="import_disk_basic">
        <sources>data/disk_basic.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="disk_basic">
            <rowmaker idmaps="*">
                <map key="rad_value" nullExpr="1e+20" />
                <map key="rad_qual" nullExpr="'?'" />
                <map key="rad_rel" nullExpr="'?'" />
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
        <column name="membership" type="integer"
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

    <data id="import_h_link">
        <sources>data/h_link.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="h_link">
             <rowmaker idmaps="*">
                 <map key="parent_object_idref" nullExpr="999999"/>
                 <map key="child_object_idref" nullExpr="999999"/>
                 <map key="membership" nullExpr="999999"/>
                 <map key="h_link_source_idref" nullExpr="999999"/>
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

    <data id="import_ident">
        <sources>data/ident.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
        <make table="ident">
             <rowmaker idmaps="*">
             </rowmaker>
        </make>
    </data>

    <table id="mes_mass_pl" onDisk="True" adql="True">
        <meta name="title">Mass measurement table</meta>
        <meta name="description">
        A list of the planetary mass measurements.

        \betawarning
        </meta>
        <primary>object_idref,mass_pl_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="mass_pl_value" type="double precision"
            ucd="phys.mass" unit="'jupiterMass'"
            tablehead="Planet Mass"
            description="Mass"
            verbLevel="1" displayHint="sf=3"/>
        <column name="mass_pl_err" type="double precision"
            ucd="stat.error;phys.mass" unit="'jupiterMass'"
            tablehead="Err. Pl. Mass"
            description="Mass error"
            verbLevel="1" displayHint="sf=3"/>
        <column name="mass_pl_qual" type="text"
            ucd="meta.code.qual;phys.mass"
            tablehead="Quality mass"
            description="Mass quality (A:best, E:worst)"
            verbLevel="1"/>
        <column name="mass_pl_rel" type="text"
            ucd="phys.mass;arith.ratio"
            tablehead="mass_rel"
            description="Mass relation defining upper / lower limit or exact
            measurement through '&lt;', '>', and '='."
            verbLevel="1"/>
        <column name="mass_pl_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="mass_source_idref"
            description="Identifier of the source of the mass
             parameter."
            required="True"
            verbLevel="1"/>
    </table>

    <data id="import_mes_mass_pl">
        <sources>data/mes_mass_pl.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_mass_pl">
             <rowmaker idmaps="*">
                 <map key="mass_pl_rel" nullExpr="'?'"/>
                 <map key="mass_pl_qual" nullExpr="'?'" />
             </rowmaker>
        </make>
    </data>
    
    <table id="mes_teff_st" onDisk="True" adql="True">
        <meta name="title">Effective temperature measurement table</meta>
        <meta name="description">
        A list of the stellar effective temperature measurements.

        \betawarning
        </meta>
        <primary>object_idref,teff_st_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="teff_st_value" type="double precision"
            ucd="phys.temperature.effective" unit="K"
            tablehead="teff_st"
            description="Object effective temperature."
            verbLevel="1" displayHint="sf=2"/>
        <column name="teff_st_err" type="double precision"
            ucd="stat.error;phys.temperature.effective" unit="K"
            tablehead="teff_st_err"
            description="Object effective temperature error."
            verbLevel="1"/>
        <column name="teff_st_qual" type="text"
            ucd="meta.code.qual;phys.temperature.effective"
            tablehead="teff_st_qual"
            description="Effective temperature quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="teff_st_source_idref" type="integer"
            ucd="meta.ref;phys.temperature.effective"
            tablehead="teff_st_source_idref"
            description="Identifier of the source of the
                effective temperature parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>

    <data id="import_mes_teff_st">
        <sources>data/mes_teff_st.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_teff_st">
             <rowmaker idmaps="*">
                 <map key="teff_st_value" nullExpr="1e+20" />
                 <map key="teff_st_err" nullExpr="1e+20" />
                 <map key="teff_st_qual" nullExpr="'?'" />
                 <map key="teff_st_source_idref" nullExpr="999999" />
             </rowmaker>
        </make>
    </data>
    
    <table id="mes_radius_st" onDisk="True" adql="True">
        <meta name="title">Radius measurement table</meta>
        <meta name="description">
        A list of the stellar radius measurements.

        \betawarning
        </meta>
        <primary>object_idref,radius_st_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the
            other tables)."
            required="True"
            verbLevel="1"/>
        <column name="radius_st_value" type="double precision"
            ucd="phys.size.radius" unit="solRad"
            tablehead="radius_st"
            description="Object radius."
            verbLevel="1" displayHint="sf=2"/>
        <column name="radius_st_err" type="double precision"
            ucd="stat.error;phys.size.radius" unit="solRad"
            tablehead="radius_st_err"
            description="Object radius error."
            verbLevel="1"/>
        <column name="radius_st_qual" type="text"
            ucd="meta.code.qual;phys.size.radius"
            tablehead="radius_st_qual"
            description="Radius quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="radius_st_source_idref" type="integer"
            ucd="meta.ref;phys.size.radius"
            tablehead="radius_st_source_idref"
            description="Identifier of the source of the
                radius parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>

    <data id="import_mes_radius_st">
        <sources>data/mes_radius_st.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_radius_st">
             <rowmaker idmaps="*">
                 <map key="radius_st_value" nullExpr="1e+20" />
                 <map key="radius_st_err" nullExpr="1e+20" />
                 <map key="radius_st_qual" nullExpr="'?'" />
                 <map key="radius_st_source_idref" nullExpr="999999" />
             </rowmaker>
        </make>
    </data>
    
    <table id="mes_mass_st" onDisk="True" adql="True">
        <meta name="title">Mass measurement table</meta>
        <meta name="description">
        A list of the stellar mass measurements.

        \betawarning
        </meta>
        <primary>object_idref,mass_st_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the
            other tables)."
            required="True"
            verbLevel="1"/>   
        <column name="mass_st_value" type="double precision"
            ucd="phys.mass" unit="solMass"
            tablehead="mass_st"
            description="Object mass."
            verbLevel="1" displayHint="sf=2"/>
        <column name="mass_st_err" type="double precision"
            ucd="stat.error;phys.mass" unit="solMass"
            tablehead="mass_st_err"
            description="Object mass error."
            verbLevel="1"/>
        <column name="mass_st_qual" type="text"
            ucd="meta.code.qual;phys.mass"
            tablehead="mass_st_qual"
            description="Mass quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="mass_st_source_idref" type="integer"
            ucd="meta.ref;phys.mass"
            tablehead="mass_st_source_idref"
            description="Identifier of the source of the
                mass parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>     
    </table>

    <data id="import_mes_mass_st">
        <sources>data/mes_mass_st.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_mass_st">
             <rowmaker idmaps="*">
                 <map key="mass_st_value" nullExpr="1e+20" />
                 <map key="mass_st_err" nullExpr="1e+20" />
                 <map key="mass_st_qual" nullExpr="'?'" />
                 <map key="mass_st_source_idref" nullExpr="999999" />
             </rowmaker>
        </make>
    </data>
    
    <table id="mes_binary" onDisk="True" adql="True">
        <meta name="title">Multiplicitz measurement table</meta>
        <meta name="description">
        A list of the stellar multiplicitz measurements.

        \betawarning
        </meta>
        <primary>object_idref,binary_flag,binary_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the
            other tables)."
            required="True"
            verbLevel="1"/>   
        <column name="binary_flag" type="text"
            ucd="meta.code.multip" 
            tablehead="binary_flag"
            description="Binary flag."
            verbLevel="1" displayHint="sf=2"/>
        <column name="binary_qual" type="text"
            ucd="meta.code.qual"
            tablehead="binary_qual"
            description="Binary quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="binary_source_idref" type="integer"
            ucd="meta.ref"
            tablehead="binary_source_idref"
            description="Identifier of the source of the
                binary_flag parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>

    <data id="import_mes_binary">
        <sources>data/mes_binary.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_binary">
             <rowmaker idmaps="*">
                 <map key="binary_flag" nullExpr="'None'" />
                 <map key="binary_qual" nullExpr="'?'" />
                 <map key="binary_source_idref" nullExpr="999999" />
             </rowmaker>
        </make>
    </data>
    
    <table id="mes_sep_phys" onDisk="True" adql="True">
        <meta name="title">Phys. separation measurement table</meta>
        <meta name="description">
        A list of the stellar phys. separation measurements.

        \betawarning
        </meta>
        <primary>object_idref,sep_phys_source_idref</primary>
        <foreignKey source="object_idref" inTable="object"
            dest="object_id" />

        <column name="object_idref" type="integer"
            ucd="meta.id;meta.main"
            tablehead="object_idref"
            description="Object key (unstable, use only for joining to the
            other tables)."
            required="True"
            verbLevel="1"/>   
        <column name="sep_phys_value" type="double precision"
            ucd="pos.angDistance" unit="arcsec"
            tablehead="Phys. separation"
            description="Phys. separation of binary."
            verbLevel="1" displayHint="sf=2"/>
        <column name="sep_phys_err" type="double precision"
            ucd="stat.error;pos.angDistance" unit="arcsec"
            tablehead="sep_phys_err"
            description="Object phys. separation error."
            verbLevel="1"/>
        <column name="sep_phys_qual" type="text"
            ucd="meta.code.qual;pos.angDistance"
            tablehead="sep_phys_qual"
            description="Phys. separation quality (A:best, E:worst)"
            verbLevel="1">
        </column>
        <column name="sep_phys_source_idref" type="integer"
            ucd="meta.ref;pos.angDistance"
            tablehead="sep_phys_source_idref"
            description="Identifier of the source of the
                sep_phys parameter."
            verbLevel="1">
              <values nullLiteral="-1"/>
        </column>
    </table>

    <data id="tables" auto="False">
        <LOOP listItems="source object star_basic planet_basic disk_basic
            h_link ident mes_mass_pl">
            <events>
                <make table="\item"/>
            </events>
        </LOOP>
        <publish sets="ivo_managed,local"/>
    </data>
    
    <data id="import_mes_sep_phys">
        <sources>data/mes_sep_phys.xml</sources>
        <!-- Data acquired using the skript data_acquisition.py. -->
        <voTableGrammar/>
           <make table="mes_sep_phys">
             <rowmaker idmaps="*">
                 <map key="sep_phys_value" nullExpr="1e+20" />
                 <map key="sep_phys_err" nullExpr="1e+20" />
                 <map key="sep_phys_qual" nullExpr="'?'" />
                 <map key="sep_phys_source_idref" nullExpr="999999" />
             </rowmaker>
        </make>
    </data>
    

    <table id="scs_summary" onDisk="True">
        <meta name="description">A view containing the most
        basic data for simple consumption by the cone search service.
        </meta>

        <LOOP>
            <csvItems>
                src, name
                object.main_id, star_name
                object.main_id, planet_name
                star_basic.coo_ra, ra
                star_basic.coo_dec, dec
                star_basic.dist_st_value, dist
                planet_basic.mass_pl_value, planet_mass
                planet_basic.mass_pl_err, planet_mass_error
            </csvItems>
            <events>
                <column original="\src" name="\name" verbLevel="5"/>
            </events>
        </LOOP>
        <column original="star_name"
            tablehead="Host Star"
            description="Name of the planet's host star"/>
        <column original="planet_name"
            tablehead="Planet"
            description="Name of the planet"
            ucd="meta.id;meta.main"/>
        <viewStatement>
            CREATE VIEW \qName AS (SELECT \colNames FROM (
            select distinct
                star_ob.main_id as star_name,
                planet_ob.main_id as planet_name,
                coo_ra as ra,
                coo_dec as dec,
                dist_st_value as dist,
                mass_pl_value as planet_mass,
                mass_pl_err as planet_mass_error
            from life_td.star_basic as s
                join life_td.h_link as slink on (parent_object_idref=s.object_idref)
                join life_td.planet_basic as p on (child_object_idref=p.object_idref)
                join life_td.object as star_ob on (s.object_idref=star_ob.object_id)
                join life_td.object as planet_ob on (p.object_idref=planet_ob.object_id))q)

        </viewStatement>
    </table>

    <data id="make_cone_view">
        <make table="scs_summary"/>
    </data>

    <service id="cone" allowed="form,scs.xml">
        <meta name="shortName">life_td cone</meta>
        <meta>
            testQuery.ra: 312.27
            testQuery.dec: 37.47
            testQuery.sr: 0.01
        </meta>
        <publish render="scs.xml" sets="ivo_managed"/>
        <publish render="form" sets="ivo_managed,local"/>

        <scsCore queriedTable="scs_summary">
            <FEED source="//scs#coreDescs"/>
            <condDesc buildFrom="planet_name"/>
            <condDesc buildFrom="planet_mass"/>
        </scsCore>
    </service>

    <service id="ex" allowed="examples">
        <meta name="title">LIFE TAP examples</meta>
        <meta name="_example" title="Filter objects by type">
            In LIFE, we have a single table for all kinds of objects
            (planets, stars, disks, multi-star systems).  They are kept in
            :taptable:`life_td.object`:

            .. tapquery::
                SELECT TOP 10 object_id, main_id FROM life_td.object
                WHERE type='st'
        </meta>
        <meta name="_example" title="All children of an object">
            Objects in LIFE are in a hierarchy (e.g., a planet belongs to a
            star).  The parent/child relationships are given in the
            :taptable:`life_td.h_link` table which you can join to all other
            tables that have an object_idref column.  For instance,
            to find (direct) children of a star, you would run:

            .. tapquery::
                SELECT main_id as Child_main_id, object_id as child_object_id
                FROM life_td.h_link
                JOIN life_td.ident as p on p.object_idref=parent_object_idref
                JOIN life_td.object on object_id=child_object_idref
                WHERE p.id = '* alf Cen'
                </meta>
        <meta name="_example" title="All parents of an object">
            Objects in LIFE are in a hierarchy (e.g., a planet belongs to a
            star).  The parent/child relationships are given in the
            :taptable:`life_td.h_link` table which you can join to all other
            tables that have an object_idref column.  For instance,
            to find (direct) parents of a star, you would run:

            .. tapquery::
                SELECT main_id as parent_main_id, object_id as parent_object_id
                FROM life_td.h_link
                JOIN life_td.ident as c on c.object_idref=child_object_idref
                JOIN life_td.object on object_id=parent_object_idref
                WHERE c.id =  '* alf Cen A'
                </meta>
        <meta name="_example" title="All specific measurements of an object">
            In LIFE, we have individual tables for all kinds of parameters
            where multiple measurements for the same object are available.
            They are kept in the tables starting with ``mes_`` e.g.
            :taptable:`life_td.mes_teff_st`:

            .. tapquery::
                SELECT *
                FROM life_td.mes_teff_st
                JOIN life_td.ident USING(object_idref)
                WHERE id = 'GJ    10'
                </meta>
        <meta name="_example" title="All basic stellar data from an object name">
            In LIFE we keep for each object the best measurements of its kind
            in the basic data table corresponding to the object type. For
            instance, to find the best measurements for the star '* alf Cen'
            you would run:

            .. tapquery::
                SELECT  *
                FROM life_td.star_basic
                JOIN life_td.ident USING(object_idref)
                WHERE id = '* alf Cen'
                </meta>
        <meta name="_example" title="All basic disk data from host name">
            In LIFE we keep for each object the best measurements of its kind
            in the basic data table corresponding to the object type. For
            instance, to find the best measurements for the disk around the
            star '* bet Cas' you would run:

            .. tapquery::
                SELECT main_id disk_main_id, object_id as disk_object_id, db.*
                FROM life_td.h_link
                JOIN life_td.disk_basic as db on db.object_idref=child_object_idref
                JOIN life_td.ident as p on p.object_idref=parent_object_idref
                JOIN life_td.object on object_id=child_object_idref
                WHERE p.id = '* bet Cas' and type='di'
                </meta>
        <meta name="_example" title="Missing reliable measurements">
            In LIFE we keep information about the quality of a measurement.
            This can serve as motivation for future observations to fill in
            knowledge gaps. For instance, to find where reliable measurements
            for the parallax are missing you would run:

            .. tapquery::
                SELECT star_ob.main_id as star_name, plx_value, plx_err,
                plx_qual, plx_source_idref
                FROM life_td.star_basic as s
                JOIN life_td.object as star_ob on
                (s.object_idref=star_ob.object_id)
                WHERE plx_value is Null or plx_qual in ('D','E') or
                plx_qual is Null
                </meta>
        <nullCore/>
    </service>

    <regSuite title="LIFE regression">
        <!-- NOTE: These tests will break per release right now because
        they use refs/idrefs; let's make them more stable when the next
        release comes around. -->

        <regTest title="LIFE form service appears to work.">
        <url parSet="form" hscs_pos="14 Her" hscs_sr="1">cone/form</url>
        <code>
            self.assertHasStrings("HD 145675c",
                "17.90", "7.679")
        </code>
        </regTest>

        <regTest title="LIFE tables appear to be in place.">
            <url parSet="TAP" QUERY="
SELECT * FROM
  life_td.planet_basic AS p
  JOIN life_td.object ON (p.object_idref=object_id)
  JOIN life_td.h_link ON (child_object_idref=object_id)
  JOIN life_td.star_basic AS s ON (parent_object_idref=s.object_idref)
WHERE
  main_id='*  14 Her b'
            ">/tap/sync</url>
            <code>
                rows = self.getVOTableRows()
                self.assertEqual(len(rows), 2)
                self.assertAlmostEqual(rows[0]["coo_ra"],
                    242.60131531625294)
                ids = set("|".join(r["ids"] for r in rows).split("|"))
                self.assertEqual(ids, {'*  14 Her b', '14 Her b',
                    'HD 145675b', 'GJ   614 b'})
            </code>
        </regTest>
    </regSuite>
</resource>
<!-- vim:et:sw=4:sta
-->
