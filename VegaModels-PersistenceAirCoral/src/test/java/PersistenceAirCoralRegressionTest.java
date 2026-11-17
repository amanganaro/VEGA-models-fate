import insilico.core.model.InsilicoModel;
import insilico.persistence_air_coral.ismPersistenceAirCoral;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the PersistenceAirCoral model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class PersistenceAirCoralRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismPersistenceAirCoral());
    }
}
