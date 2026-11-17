import insilico.core.model.InsilicoModel;
import insilico.readybio_irfmn.ismReadyBioIRFMN;
import regression.ModelRegressionTest;

import java.util.List;

/**
 * Regression test of the ReadyBioIRFMN model, see {@link ModelRegressionTest}.
 * Run with {@code mvn test -Pregression}; regenerate the references with {@code -Dgolden.update=true}.
 */
public class ReadyBioIRFMNRegressionTest extends ModelRegressionTest {

    @Override
    protected List<InsilicoModel> getModels() throws Exception {
        return List.of(new ismReadyBioIRFMN());
    }
}
