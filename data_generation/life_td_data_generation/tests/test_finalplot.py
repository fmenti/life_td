from astropy.table import Table
from utils.analysis.finalplot import starcat_distribution_plot

cat1 = Table(
    [
        ["M", "M","G","K","F","A","M", "M","G","K","F","A"],
        [10,19,5,3,2,15,2,17,5,14,2,1],
    ],
    names=["class_temp", "dist_st_value"],
    dtype=[object, float],
)
cat2 = Table(
    [
        [ "M","G","K","F","G","M", "M","G","K","M","A"],
        [5,3,2,15,2,17,5,14,2,1,15],
    ],
    names=["class_temp", "dist_st_value"],
    dtype=[object, float],
)
starcat_distribution_plot([cat1,cat2],['test_data1','test_data2'])
