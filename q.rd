<resource schema="life">
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
		<primary>main_id</primary>
		<column name="main_id" type="text"
			ucd="meta.id;meta.main"
			tablehead="main_id"
			description="Main stellar identifier."
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
	</table>    
    
	<data id="import_main_id">
		<sources pattern="data/simbad_main_id.csv"/>
			<!-- aquired via TopCat TAP SIMBAD query:
			SELECT TOP 10 main_id,ra,dec
			FROM basic
			WHERE otype='*' -->
		<csvGrammar/>
		<make table="star_basic">
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
