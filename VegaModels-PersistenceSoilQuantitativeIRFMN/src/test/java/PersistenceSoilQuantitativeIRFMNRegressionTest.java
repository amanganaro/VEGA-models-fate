import insilico.core.model.InsilicoModel;
import insilico.persistence_soil_quantitative_irfmn.ismPersistenceSoilQuantitativeIrfmn;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the PersistenceSoilQuantitativeIRFMN model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class PersistenceSoilQuantitativeIRFMNRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismPersistenceSoilQuantitativeIrfmn());
    }
}
