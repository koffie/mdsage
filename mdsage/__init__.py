# Add the import for which you want to give a direct access

try:
    from .ultimate_question import answer_to_ultimate_question
except ImportError:
    pass

from .cuspidal_classgroup import (
    intersection,
    cuspidal_rational_subgroup_mod_rational_cuspidal_subgroup,
    upper_bound_index_cusps_in_JG_torsion,
)

from .kamiennys_criterion import KamiennyCriterion, matrix_modp

from .maartens_sage_functions import (
    congruence_groups_between_gamma0_and_gamma1,
    count_points_J_H,
    counts,
    gonality_lower_bound,
    cuspidal_integral_structure_matrix,
    generators_of_subgroups_of_unit_group,
    diamond_operator,
    positive_part,
    rational_cuspidal_classgroup,
    tate_normal_form,
)

from .modular_unit_divisors import cuspsums_on_Gamma1_of_degree

from .quadratic_class_numbers import (
    class_number,
    class_numbers,
    cm_orders2,
    small_class_number_fundamental_discriminants,
    small_class_number_discriminants,
)

from .ramification import (
    small_ramification,
    atkin_lehner_ramification_degree,
    atkin_lehner_divisors,
)

from .congruence_subgroups import intermediate_modular_groups, congruence_subgroup_repr

from .canonical_rings import subhyperelliptic_primes

from .modular_degrees_oldforms import (
    product_isogeny_map,
    modular_symbol_elliptic_curves,
    modular_symbol_elliptic_curves_range,
    modular_symbol_elliptic_curves_divisors,
    degree_pairing,
    degree_quadratic_form,
)
