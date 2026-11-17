import insilico.core.model.InsilicoModel;
import insilico.hydrolysis_coral.ismHydrolysisCoral;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the Hydrolysis (IRFMN-CORAL) model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class HydrolysisCoralRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismHydrolysisCoral());
    }
}
