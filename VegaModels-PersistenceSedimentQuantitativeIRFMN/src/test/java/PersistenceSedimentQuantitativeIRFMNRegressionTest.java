import insilico.core.model.InsilicoModel;
import insilico.persistence_sediment_quantitative_irfmn.ismPersistenceSedimentQuantitativeIrfmn;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the PersistenceSedimentQuantitativeIRFMN model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class PersistenceSedimentQuantitativeIRFMNRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismPersistenceSedimentQuantitativeIrfmn());
    }
}
