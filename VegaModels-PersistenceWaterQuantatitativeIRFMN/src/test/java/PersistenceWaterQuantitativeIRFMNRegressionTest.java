import insilico.core.model.InsilicoModel;
import insilico.persistence_quantative_water_irfmn.ismPersistenceWaterQuantitativeIrfmn;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the PersistenceWaterQuantatitativeIRFMN model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class PersistenceWaterQuantitativeIRFMNRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismPersistenceWaterQuantitativeIrfmn());
    }
}
