from typing import NamedTuple

from constants import FRUIT_DEVELOPMENT_STAGES_NUM


class TomatoModel(NamedTuple):
    carbohydrate_amount_Buf: float
    carbohydrate_amount_Fruits: [float]*FRUIT_DEVELOPMENT_STAGES_NUM
    number_Fruits: [float]*FRUIT_DEVELOPMENT_STAGES_NUM
    carbohydrate_amount_Leaf: float
    carbohydrate_amount_Stem: float
    dry_matter_Har: float
    sum_canopy_t: float
    last_24_canopy_t: float

