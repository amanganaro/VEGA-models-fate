import insilico.core.model.InsilicoModel;
import insilico.persistence_water_irfmn.ismPersistenceWaterIrfmn;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the PersistenceWaterIRFMN model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class PersistenceWaterIRFMNRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismPersistenceWaterIrfmn());
    }
}
