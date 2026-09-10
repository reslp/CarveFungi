'''
  Modified version of the CarveMe algorithm. Code extracted from: https://github.com/cdanielmachado/carveme								
  D. Machado et al, "Fast automated reconstruction of genome-scale metabolic models for microbial species and communities", Nucleic Acids Research, gky537, 2018..
  DOI: https://doi.org/10.1093/nar/gky537							
  Developed at the European Molecular Biology Laboratory (2017-2018). Released under an Apache License. 
   
'''
import cobra
from pyscipopt import Model as SCIPModel, quicksum
import operator
import warnings


class solution:
    obj=None
    sol=None
    def __init__(self, obj, sol):
        self.obj = obj
        self.sol = sol
    def get_variable_names(self):
        return self.sol.keys()
    def get_solution(self):
        return self.sol
    def get_obj(self):
        return self.obj

def generate_soln_pool(solver, relative_gap=0.001):
    """Return stored SCIP solutions within 0.1% of the best incumbent score.

    Unlike CPLEX's populate_solution_pool, this does not actively enumerate
    alternatives. SCIP may retain only one qualifying solution.
    """
    solver.optimize()
    status = str(solver.getStatus())
    candidates = sorted(solver.getSols(), key=solver.getSolObjVal, reverse=True)
    if not candidates:
        raise RuntimeError(f"SCIP returned no feasible reconstruction (status: {status}).")
    if status != "optimal":
        warnings.warn(
            f"SCIP status: {status}; using feasible incumbents, not proven optimal "
            f"(MIP gap: {solver.getGap():.6g}).",
            RuntimeWarning,
        )

    best = solver.getSolObjVal(candidates[0])
    tolerance = relative_gap * abs(best) + 1e-9
    variables = solver.getVars(transformed=False)
    sol_pool = []
    for candidate in candidates:
        objective = solver.getSolObjVal(candidate)
        if best - objective <= tolerance:
            values = {var.name: solver.getSolVal(candidate, var) for var in variables}
            sol_pool.append(solution(objective, values))

    print(f"SCIP status: {status}; retained {len(sol_pool)} solutions; best score: {best:.10g}")
    return sol_pool


def minmax_reduction(model, scores, min_growth=0.1, min_atpm=0.1, eps=1e-5, bigM=1e3, default_score=-1.0,
                     uptake_score=1, ref_reactions=[], ref_score=0.0, feast=1e-6, opti=1e-5):
    """ Apply minmax reduction algorithm (MILP).
    Computes a binary reaction vector that optimizes the agreement with reaction scores (maximizes positive scores,
    and minimizes negative scores). It generates a fully connected reaction network (i.e. all reactions must be able
    to carry some flux).
    Args:
        model (CBModel): universal model
        scores (dict): reaction scores
        min_growth (float): minimal growth constraint
        min_atpm (float): minimal maintenance ATP constraint
        eps (float): minimal flux required to consider leaving the reaction in the model
        bigM (float): maximal reaction flux
        default_score (float): penalty score for reactions without an annotation score (default: -1.0).
        uptake_score (float): penalty score for using uptake reactions (default: 0.0).
    Returns:
        Solution: optimization result
    """
    problem = SCIPModel("CarveFungi")
    problem.setRealParam("numerics/dualfeastol", opti)
    problem.setRealParam("numerics/feastol", feast)
    problem.setMaximize()

    variables = {}
    for reaction in model.reactions:
        variables[reaction.id] = problem.addVar(
            name=reaction.id, lb=reaction.lower_bound, ub=reaction.upper_bound)
    for metabolite in model.metabolites:
        balance = quicksum(
            reaction.metabolites[metabolite] * variables[reaction.id]
            for reaction in sorted(metabolite.reactions, key=lambda r: r.id)
        )
        problem.addCons(balance == 0, name="mass_" + metabolite.id)

    scores = scores.copy()

    reactions = []
    for r_id in model.reactions:
        if r_id.id not in reactions and not r_id.id.endswith('_E'):
            reactions.append(r_id.id)


    # R_UF01847_CE is the ATP maintenance reaction
    # if the default score is lower than 0 set all the reactions to the default score except for the exchange reactions
    # and the ATP maintenance reaction
    if default_score != 0:
        for r_id in model.reactions:
            if r_id.id not in ref_reactions and not r_id.id.endswith('_E') and r_id.id != 'UF01847_CE' and r_id.id not in reactions:
                scores[r_id.id] = default_score
                reactions.append(r_id.id)

    if ref_score != 0:
        for r_id in ref_reactions and r_id not in reactions:
            if r_id != 'UF01847_CE':
                scores[r_id] = ref_score
                reactions.append(r_id)

    problem.addCons(variables['BIOMASS'] >= min_growth, name='min_growth')
    problem.addCons(variables['UF01847_CE'] >= min_atpm, name='min_atpm')

    neg_vars = []
    pos_vars = []

    for r_id in reactions:
        rxn = model.reactions.get_by_id(r_id)
        if rxn.lower_bound < 0:
            y_r = 'yr_' + r_id
            variables[y_r] = problem.addVar(name=y_r, vtype="B", obj=float(scores[r_id]))
            neg_vars.append(y_r)
        if rxn.upper_bound > 0:
            y_f = 'yf_' + r_id
            variables[y_f] = problem.addVar(name=y_f, vtype="B", obj=float(scores[r_id]))
            pos_vars.append(y_f)

    if uptake_score != 0:
        for r_id in model.reactions:
            if r_id.id.endswith('_E'):
                name = 'y_' + r_id.id
                variables[name] = problem.addVar(name=name, vtype="B", obj=uptake_score)

    for r_id in reactions:
        y_r, y_f = 'yr_' + r_id, 'yf_' + r_id
        flux = variables[r_id]
        if y_r in neg_vars and y_f in pos_vars:
            forward, reverse = variables[y_f], variables[y_r]
            problem.addCons(flux - eps * forward + bigM * reverse >= 0, name='lb_' + r_id)
            problem.addCons(flux - bigM * forward + eps * reverse <= 0, name='ub_' + r_id)
            problem.addCons(forward + reverse <= 1, name='rev_' + r_id)
        elif y_f in pos_vars:
            forward = variables[y_f]
            problem.addCons(flux - eps * forward >= 0, name='lb_' + r_id)
            problem.addCons(flux - bigM * forward <= 0, name='ub_' + r_id)
        elif y_r in neg_vars:
            reverse = variables[y_r]
            problem.addCons(flux + bigM * reverse >= 0, name='lb_' + r_id)
            problem.addCons(flux + eps * reverse <= 0, name='ub_' + r_id)

    if uptake_score != 0:
        for r_id in model.reactions:
            if r_id.id.endswith('_E'):
                problem.addCons(variables[r_id.id] + bigM * variables['y_' + r_id.id] >= 0,
                                name='lb_' + r_id.id)
    solutions = generate_soln_pool(problem)

    return solutions


def carve_model(model, reaction_scores, min_growth=0.1, min_atpm=0.1,
                default_score=-0.1, uptake_score=1.0, eps=1e-5, feast=1e-6, opti=1e-6):
    """ Reconstruct a metabolic model using the CarveMe approach.
    Args:
        model (CBModel): universal model
        reaction_scores (pandas.DataFrame): reaction scores
        outputfile (str): write model to SBML file (optional)
        flavor (str): SBML flavor ('cobra' or 'fbc2', optional)
        inplace (bool): Change model in place (default: True)
        default_score (float): penalty for non-annotated intracellular reactions (default: -1.0)
        uptake_score (float): penalty for utilization of extracellular compounds (default: 0.0)
    Returns:
        CBModel: reconstructed model
    """

    scores = reaction_scores

    if default_score is None:
        default_score = -1.0

    if uptake_score is None:
        uptake_score = 0.0


    reconstructed_model = model.copy()

    solutions = minmax_reduction(reconstructed_model, scores, eps=eps, default_score=default_score, uptake_score=uptake_score,
                                 min_growth=min_growth, min_atpm=min_atpm, feast=feast, opti=opti)

    solutions = select_best_models(solutions)

    models = []
    i = 0
    objective = []
    try:
        for sol, obj in solutions.items():
            reconstructed_model = model.copy()
            print("Solution " + str(i))
            i = i + 1
            inactive = inactive_reactions(reconstructed_model, sol)
            reconstructed_model.remove_reactions(inactive)
            reconstructed_model, removed = cobra.manipulation.delete.prune_unused_metabolites(reconstructed_model)
            reconstructed_model.objective = "BIOMASS"
            print(reconstructed_model.optimize())
            models.append(reconstructed_model)
            objective.append(obj)
    except:
        raise

    return objective, models


def select_best_models(solutions):
    best_candidates = {}
    for sol in solutions:
        best_candidates[sol] = sol.get_obj()
    best_solutions = dict(sorted(best_candidates.items(), key=operator.itemgetter(1), reverse=True)[0:5])
    return best_solutions


def inactive_reactions(model, solution):
    inactive = []

    internal = [r_id.id for r_id in model.reactions if not r_id.id.endswith('_E')]
    external = [r_id.id for r_id in model.reactions if r_id.id.endswith('_E')]

    remove=[]

    solu_values=solution.get_solution()
    for r_id in internal:
        threshold = 0.5
        try:
            if "yf_" + r_id in solu_values and solu_values['yf_' + r_id] < threshold and "yr_" + r_id in solu_values and \
                    solu_values['yr_' + r_id] < threshold:
                inactive.append(r_id)
            if "yf_" + r_id not in solu_values and "yr_" + r_id in solu_values and solu_values[
                'yr_' + r_id] < threshold:
                inactive.append(r_id)
            if "yf_" + r_id in solu_values and solu_values[
                'yf_' + r_id] < threshold and "yr_" + r_id not in solu_values:
                inactive.append(r_id)
        except:
            raise

    inactive_ext = []
    for r_id in external:
        try:
            neighbors=[]
            for metabolite in model.reactions.get_by_id(r_id).metabolites:
                for rxn in metabolite.reactions:
                    neighbors.append(rxn.id)
            if len(set(neighbors) - set(inactive)) == 1:
                inactive_ext.append(r_id)
        except:
            print(r_id, m_id, "is there something wrong?")
    return inactive + inactive_ext+remove
